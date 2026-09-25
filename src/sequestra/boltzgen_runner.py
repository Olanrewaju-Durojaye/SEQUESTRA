"""Prepare and execute the resumable BoltzGen generation stage."""

from __future__ import annotations

import json
import os
from pathlib import Path
import signal
import subprocess
import threading
from typing import Any, Callable

from .run_state import (
    finish_stage, interrupt_stage, record_partial_stage_success,
    recover_interrupted, start_stage,
)


def _load_json_yaml(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def prepare_design_specification(project_dir: Path) -> Path:
    root = project_dir.expanduser().resolve()
    specification = _load_json_yaml(root / "configuration/project_specification.yaml")
    design = specification["design"]
    ligand = specification["ligand"]
    ccd = str(ligand["ccd"])
    minimum = int(design["minimum_length"])
    maximum = int(design["maximum_length"])
    path = root / "boltzgen/design_specification.yaml"
    content = (
        "entities:\n"
        "  - protein:\n"
        "      id: A\n"
        f"      sequence: {minimum}..{maximum}\n"
        "  - ligand:\n"
        "      id: B\n"
        f"      ccd: {ccd}\n"
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def resolve_boltzgen_executable(config_path: Path) -> Path:
    config = json.loads(config_path.expanduser().read_text(encoding="utf-8"))
    environment = Path(config["environments"]["boltzgen"]).expanduser()
    executable = environment / "bin/boltzgen"
    if not executable.is_file():
        raise FileNotFoundError(f"BoltzGen executable not found: {executable}")
    return executable.resolve()


def _resource_settings(config_path: Path) -> tuple[int, float, float]:
    """Return conservative BoltzGen worker and host-memory safeguards."""
    config = json.loads(config_path.expanduser().read_text(encoding="utf-8"))
    settings = config.get("boltzgen", {})
    workers = int(settings.get("num_workers", 0))
    minimum_gib = float(settings.get("minimum_available_ram_gib", 6.0))
    poll_seconds = float(settings.get("memory_poll_seconds", 2.0))
    if workers < 0:
        raise ValueError("boltzgen.num_workers must be zero or greater")
    if minimum_gib < 1:
        raise ValueError("boltzgen.minimum_available_ram_gib must be at least 1")
    if poll_seconds <= 0:
        raise ValueError("boltzgen.memory_poll_seconds must be greater than zero")
    return workers, minimum_gib, poll_seconds


def _available_memory_bytes() -> int | None:
    """Read Linux MemAvailable without adding a runtime dependency."""
    try:
        for line in Path("/proc/meminfo").read_text(encoding="utf-8").splitlines():
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        return None
    return None


def _terminate_process_group(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except (ProcessLookupError, PermissionError):
        process.terminate()


def build_boltzgen_command(
    project_dir: Path,
    *,
    config_path: Path,
    smoke_test: bool,
    devices: int = 1,
) -> list[str]:
    root = project_dir.expanduser().resolve()
    specification = _load_json_yaml(root / "configuration/project_specification.yaml")
    requested = int(specification["design"]["number_of_designs"])
    number = min(2, requested) if smoke_test else requested
    design_spec = prepare_design_specification(root)
    executable = resolve_boltzgen_executable(config_path)
    workers, _, _ = _resource_settings(config_path)
    return [
        str(executable), "run", str(design_spec),
        "--output", str(root / "boltzgen/run"),
        "--protocol", "protein-small_molecule",
        "--num_designs", str(number),
        "--budget", str(number),
        "--config", "analysis", "num_processes=1",
        "--devices", str(devices),
        "--num_workers", str(workers),
        "--reuse",
    ]


def _expected_outputs(project_dir: Path) -> list[Path]:
    output = project_dir.expanduser().resolve() / "boltzgen/run"
    return [
        output / "final_ranked_designs/all_designs_metrics.csv",
        output / "intermediate_designs_inverse_folded/aggregate_metrics_analyze.csv",
    ]


def run_boltzgen(
    project_dir: Path,
    *,
    config_path: Path,
    smoke_test: bool = False,
    devices: int = 1,
    emit: Callable[[str], None] = print,
) -> int:
    root = project_dir.expanduser().resolve()
    command = build_boltzgen_command(
        root, config_path=config_path, smoke_test=smoke_test, devices=devices
    )
    recover_interrupted(root)
    workers, minimum_gib, poll_seconds = _resource_settings(config_path)
    available = _available_memory_bytes()
    floor_bytes = int(minimum_gib * 1024**3)
    if available is not None and available < floor_bytes:
        emit(
            f"BoltzGen not started: {available / 1024**3:.1f} GiB RAM is available; "
            f"the configured safety floor is {minimum_gib:.1f} GiB."
        )
        return 75
    log_path = root / "logs/boltzgen.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start_stage(root, "boltzgen", command=command)
    emit("Executing: " + " ".join(command))
    emit(f"Log: {log_path}")
    emit(
        f"Resource safeguards: DataLoader workers={workers}; "
        f"minimum available RAM={minimum_gib:.1f} GiB."
    )
    process: subprocess.Popen[str] | None = None
    stop_monitor = threading.Event()
    memory_stop = threading.Event()
    monitor: threading.Thread | None = None
    try:
        with log_path.open("a", encoding="utf-8") as log:
            process = subprocess.Popen(
                command,
                cwd=root,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                start_new_session=True,
            )
            def monitor_memory() -> None:
                consecutive_low = 0
                while not stop_monitor.wait(poll_seconds):
                    current = _available_memory_bytes()
                    if current is None:
                        continue
                    consecutive_low = consecutive_low + 1 if current < floor_bytes else 0
                    if consecutive_low >= 3:
                        memory_stop.set()
                        _terminate_process_group(process)
                        return

            monitor = threading.Thread(target=monitor_memory, daemon=True)
            monitor.start()
            assert process.stdout is not None
            for line in process.stdout:
                log.write(line)
                log.flush()
                emit(line.rstrip())
            process.stdout.close()
            return_code = process.wait()
    except KeyboardInterrupt:
        if process is not None:
            _terminate_process_group(process)
            process.wait()
        interrupt_stage(root, "boltzgen", message="interrupted by user")
        emit("BoltzGen interrupted. Rerun the same command to continue with --reuse.")
        return 130
    finally:
        stop_monitor.set()
        if monitor is not None:
            monitor.join(timeout=max(1.0, poll_seconds * 2))
    if memory_stop.is_set():
        message = (
            f"BoltzGen stopped safely because available RAM remained below "
            f"{minimum_gib:.1f} GiB; completed artifacts were preserved for --reuse"
        )
        finish_stage(root, "boltzgen", succeeded=False, message=message)
        emit(message)
        emit("Close other applications or increase swap, then run 'sequestra run' again.")
        return 75
    if return_code != 0:
        finish_stage(
            root, "boltzgen", succeeded=False,
            message=f"BoltzGen exited with status {return_code}",
        )
        return return_code
    existing = [path for path in _expected_outputs(root) if path.is_file()]
    if not existing:
        finish_stage(
            root, "boltzgen", succeeded=False,
            message="BoltzGen returned success but expected metrics outputs were not found",
        )
        emit("Output validation failed; see the BoltzGen log.")
        return 2
    relative_outputs = [str(path.relative_to(root)) for path in existing]
    if smoke_test:
        record_partial_stage_success(
            root, "boltzgen",
            message="two-design smoke test completed; full requested campaign remains pending",
            outputs=relative_outputs,
        )
        emit("Smoke test completed. The full BoltzGen stage remains pending.")
    else:
        finish_stage(
            root, "boltzgen", succeeded=True,
            message="requested BoltzGen campaign completed",
            outputs=relative_outputs,
        )
        emit("BoltzGen stage completed and checkpointed.")
    return 0
