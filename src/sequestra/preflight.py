"""Read-only preflight checks for the SEQUESTRA workstation integration."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any


GIB = 1024**3
SCIENTIFIC_PACKAGES = {
    "boltzgen": "boltzgen",
    "boltz2": "boltz",
    "catpred": "catpred",
}


def _result(name: str, status: str, detail: str) -> dict[str, str]:
    return {"name": name, "status": status, "detail": detail}


def _load_config(path: Path) -> dict[str, Any]:
    resolved = path.expanduser().resolve()
    with resolved.open(encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("schema_version") != 1:
        raise ValueError("unsupported or missing configuration schema_version")
    if not isinstance(config.get("environments"), dict):
        raise ValueError("configuration lacks an environments object")
    return config


def _package_version(python: Path, distribution: str) -> tuple[str, str]:
    """Read package metadata in a target environment without importing the tool."""
    code = (
        "import importlib.metadata as m; "
        f"print(m.version({distribution!r}))"
    )
    completed = subprocess.run(
        [str(python), "-c", code],
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    if completed.returncode == 0:
        return "pass", completed.stdout.strip()
    return "warn", "distribution metadata unavailable"


def _memory_available_gib() -> float | None:
    meminfo = Path("/proc/meminfo")
    if not meminfo.is_file():
        return None
    for line in meminfo.read_text(encoding="utf-8").splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) * 1024 / GIB
    return None


def _gpu_info() -> tuple[str, int | None, str]:
    executable = shutil.which("nvidia-smi")
    if not executable:
        return "warn", None, "nvidia-smi was not found"
    completed = subprocess.run(
        [
            executable,
            "--query-gpu=name,memory.total,driver_version",
            "--format=csv,noheader,nounits",
        ],
        capture_output=True,
        text=True,
        timeout=20,
        check=False,
    )
    if completed.returncode != 0 or not completed.stdout.strip():
        return "warn", None, "GPU query failed"
    first = completed.stdout.strip().splitlines()[0]
    parts = [part.strip() for part in first.split(",")]
    try:
        memory_mib = int(parts[1])
    except (IndexError, ValueError):
        memory_mib = None
    return "pass", memory_mib, first


def run_preflight(config_path: Path) -> dict[str, Any]:
    """Run metadata and resource checks only; never start a scientific model."""
    config = _load_config(config_path)
    checks: list[dict[str, str]] = []

    conda = Path(config["conda_executable"]).expanduser()
    checks.append(
        _result(
            "conda executable",
            "pass" if conda.is_file() and os.access(conda, os.X_OK) else "fail",
            str(conda),
        )
    )

    environments: dict[str, str] = config["environments"]
    for name, raw_path in environments.items():
        root = Path(raw_path).expanduser()
        python = root / "bin/python"
        exists = root.is_dir() and python.is_file() and os.access(python, os.X_OK)
        checks.append(
            _result(f"environment {name}", "pass" if exists else "fail", str(root))
        )
        distribution = SCIENTIFIC_PACKAGES.get(name)
        if exists and distribution:
            status, version = _package_version(python, distribution)
            checks.append(_result(f"package {distribution}", status, version))

    for name, raw_path in config.get("repositories", {}).items():
        path = Path(raw_path).expanduser()
        checks.append(
            _result(f"repository {name}", "pass" if path.is_dir() else "fail", str(path))
        )

    for name in ("dlkcat", "catpred"):
        adapter = config.get("predictors", {}).get(name, {})
        repository = Path(adapter.get("repository", "")).expanduser()
        if adapter.get("adapter") == "native":
            required = ([repository / "DeeplearningApproach/Code/example/prediction_for_input.py", repository / "DeeplearningApproach/Data", repository / "DeeplearningApproach/Results"] if name == "dlkcat" else [repository / "predict.py", repository / "scripts/create_pdbrecords.py", Path(adapter.get("data_root", "")).expanduser()])
            missing = [str(path) for path in required if not path.exists()]
            if name == "catpred" and not adapter.get("checkpoint_dir") and required[-1].is_dir():
                found = any(path.is_dir() and path.name.lower() == "km" and any(item.is_file() for item in path.rglob("*")) for path in required[-1].rglob("km"))
                if not found: missing.append(f"Km checkpoints under {required[-1]}")
            elif name == "catpred" and adapter.get("checkpoint_dir") and not Path(adapter["checkpoint_dir"]).expanduser().is_dir():
                missing.append(str(Path(adapter["checkpoint_dir"]).expanduser()))
            valid = repository.is_dir() and not missing
            detail = str(repository) if valid else "missing: " + ", ".join(missing)
        else:
            command = adapter.get("command"); output = adapter.get("output_csv")
            valid = repository.is_dir() and isinstance(command, list) and bool(command) and bool(output)
            detail = str(repository)
        checks.append(_result(f"direct {name} adapter", "pass" if valid else "fail", detail))

    boltz_executable = Path(environments["boltz2"]).expanduser() / "bin/boltz"
    checks.append(_result("direct Boltz-2 executable", "pass" if boltz_executable.is_file() and os.access(boltz_executable, os.X_OK) else "fail", str(boltz_executable)))

    minima = config.get("minimum_resources", {})
    repository = Path(config["repositories"]["sequestra"]).expanduser()
    disk_target = repository if repository.exists() else repository.parent
    free_disk = shutil.disk_usage(disk_target).free / GIB
    disk_minimum = float(minima.get("free_disk_gib", 0))
    checks.append(
        _result(
            "free disk",
            "pass" if free_disk >= disk_minimum else "warn",
            f"{free_disk:.1f} GiB available; minimum {disk_minimum:.1f} GiB",
        )
    )

    available_ram = _memory_available_gib()
    ram_minimum = float(minima.get("available_ram_gib", 0))
    checks.append(
        _result(
            "available RAM",
            "pass" if available_ram is not None and available_ram >= ram_minimum else "warn",
            "unavailable" if available_ram is None else f"{available_ram:.1f} GiB available; minimum {ram_minimum:.1f} GiB",
        )
    )

    gpu_status, gpu_memory, gpu_detail = _gpu_info()
    gpu_minimum = int(minima.get("gpu_memory_mib", 0))
    if gpu_memory is not None and gpu_memory < gpu_minimum:
        gpu_status = "warn"
    checks.append(_result("GPU", gpu_status, gpu_detail))

    active_prefix = os.environ.get("CONDA_PREFIX", "")
    expected_prefix = str(Path(environments["sequestra"]).expanduser())
    checks.append(
        _result(
            "active SEQUESTRA environment",
            "pass" if active_prefix == expected_prefix else "warn",
            active_prefix or "no active Conda environment",
        )
    )

    statuses = {item["status"] for item in checks}
    overall = "fail" if "fail" in statuses else "warn" if "warn" in statuses else "pass"
    return {
        "report_schema_version": 1,
        "read_only": True,
        "scientific_models_executed": False,
        "overall_status": overall,
        "checks": checks,
    }


def print_report(report: dict[str, Any]) -> None:
    print("SEQUESTRA preflight audit")
    print("Scientific models executed: no")
    for item in report["checks"]:
        print(f"[{item['status'].upper():4}] {item['name']}: {item['detail']}")
    print(f"Overall status: {report['overall_status'].upper()}")


def write_json_report(report: dict[str, Any], path: Path) -> None:
    destination = path.expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
