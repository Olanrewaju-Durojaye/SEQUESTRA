"""Atomic, resumable workflow state for SEQUESTRA projects."""

from __future__ import annotations

from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile
from typing import Any

from .decision_rules import CATALYTIC_REFERENCE

RUN_STATE_RELATIVE_PATH = Path("configuration/run_state.json")
STAGE_ORDER = (
    "boltzgen",
    "shortlist",
    "catalytic_liability",
    "boltz2_affinity",
    "structural_confidence",
    "reports",
)
VALID_STATUSES = {"pending", "running", "completed", "failed", "interrupted", "skipped"}


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _state_path(project_dir: Path) -> Path:
    return project_dir.expanduser().resolve() / RUN_STATE_RELATIVE_PATH


def _atomic_json_write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=path.parent
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
        directory_descriptor = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory_descriptor)
        finally:
            os.close(directory_descriptor)
    finally:
        if temporary.exists():
            temporary.unlink()


def initialize_run_state(
    project_dir: Path,
    *,
    reference_mode: str,
    catalytic_liability_enabled: bool,
    overwrite: bool = False,
) -> dict[str, Any]:
    path = _state_path(project_dir)
    if path.exists() and not overwrite:
        raise FileExistsError(f"run state already exists: {path}")
    liability_active = reference_mode == CATALYTIC_REFERENCE or catalytic_liability_enabled
    created = _now()
    stages = []
    for name in STAGE_ORDER:
        skipped = name == "catalytic_liability" and not liability_active
        stages.append({
            "name": name,
            "status": "skipped" if skipped else "pending",
            "attempts": [],
            "completed_at": created if skipped else None,
            "skip_reason": "disabled_for_non_catalytic_reference" if skipped else None,
        })
    state = {
        "schema_version": 1,
        "created_at": created,
        "updated_at": created,
        "reference_mode": reference_mode,
        "resume_policy": {
            "completed_stages_are_immutable": True,
            "resume_from_first_incomplete_stage": True,
            "running_stage_after_restart_becomes_interrupted": True,
        },
        "stages": stages,
    }
    _atomic_json_write(path, state)
    return state


def initialize_from_project(project_dir: Path) -> dict[str, Any]:
    root = project_dir.expanduser().resolve()
    specification_path = root / "configuration/project_specification.yaml"
    if not specification_path.is_file():
        raise FileNotFoundError(f"project specification not found: {specification_path}")
    specification = json.loads(specification_path.read_text(encoding="utf-8"))
    mode = specification["reference"]["mode"]
    enabled = bool(specification["decision_rules"]["catalytic_liability_enabled"])
    return initialize_run_state(
        root, reference_mode=mode, catalytic_liability_enabled=enabled
    )


def load_run_state(project_dir: Path) -> dict[str, Any]:
    path = _state_path(project_dir)
    if not path.is_file():
        raise FileNotFoundError(
            f"run state not found: {path}; run 'sequestra workflow-init --project PROJECT'"
        )
    try:
        state = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError(f"run state is not valid JSON: {path}") from error
    if state.get("schema_version") != 1:
        raise ValueError("unsupported run-state schema")
    stages = state.get("stages")
    if not isinstance(stages, list) or [row.get("name") for row in stages] != list(STAGE_ORDER):
        raise ValueError("run state has an invalid stage sequence")
    if any(row.get("status") not in VALID_STATUSES for row in stages):
        raise ValueError("run state contains an invalid stage status")
    return state


def save_run_state(project_dir: Path, state: dict[str, Any]) -> None:
    state["updated_at"] = _now()
    _atomic_json_write(_state_path(project_dir), state)


def recover_interrupted(project_dir: Path) -> dict[str, Any]:
    """Convert stale running work to interrupted so it can be resumed safely."""
    state = load_run_state(project_dir)
    changed = False
    for stage in state["stages"]:
        if stage["status"] == "running":
            stage["status"] = "interrupted"
            attempts = stage.get("attempts", [])
            if attempts and attempts[-1].get("finished_at") is None:
                attempts[-1]["finished_at"] = _now()
                attempts[-1]["outcome"] = "interrupted"
                attempts[-1]["message"] = "process ended before stage completion"
            changed = True
    if changed:
        save_run_state(project_dir, state)
    return state


def next_resumable_stage(state: dict[str, Any]) -> str | None:
    for stage in state["stages"]:
        if stage["status"] not in {"completed", "skipped"}:
            return str(stage["name"])
    return None


def start_stage(project_dir: Path, stage_name: str, *, command: list[str] | None = None) -> dict[str, Any]:
    state = recover_interrupted(project_dir)
    expected = next_resumable_stage(state)
    if expected is None:
        raise ValueError("workflow is already complete")
    if stage_name != expected:
        raise ValueError(f"cannot start {stage_name}; next resumable stage is {expected}")
    stage = next(row for row in state["stages"] if row["name"] == stage_name)
    if stage["status"] == "completed":
        raise ValueError(f"completed stage cannot be rerun: {stage_name}")
    started = _now()
    stage["status"] = "running"
    stage["attempts"].append({
        "attempt": len(stage["attempts"]) + 1,
        "started_at": started,
        "finished_at": None,
        "outcome": None,
        "command": command,
        "message": None,
    })
    save_run_state(project_dir, state)
    return state


def finish_stage(
    project_dir: Path,
    stage_name: str,
    *,
    succeeded: bool,
    message: str | None = None,
    outputs: list[str] | None = None,
) -> dict[str, Any]:
    state = load_run_state(project_dir)
    stage = next((row for row in state["stages"] if row["name"] == stage_name), None)
    if stage is None:
        raise ValueError(f"unknown stage: {stage_name}")
    if stage["status"] != "running" or not stage["attempts"]:
        raise ValueError(f"stage is not running: {stage_name}")
    finished = _now()
    stage["status"] = "completed" if succeeded else "failed"
    stage["attempts"][-1].update({
        "finished_at": finished,
        "outcome": stage["status"],
        "message": message,
        "outputs": outputs or [],
    })
    if succeeded:
        stage["completed_at"] = finished
    save_run_state(project_dir, state)
    return state


def skip_stage(
    project_dir: Path,
    stage_name: str,
    *,
    reason: str,
    outputs: list[str] | None = None,
) -> dict[str, Any]:
    """Mark the next pending stage as deliberately skipped with an auditable reason."""
    state = recover_interrupted(project_dir)
    expected = next_resumable_stage(state)
    if expected != stage_name:
        raise ValueError(
            f"cannot skip {stage_name}; next resumable stage is {expected or 'none'}"
        )
    stage = next(row for row in state["stages"] if row["name"] == stage_name)
    finished = _now()
    stage.update({
        "status": "skipped",
        "completed_at": finished,
        "skip_reason": reason,
        "skip_outputs": outputs or [],
    })
    save_run_state(project_dir, state)
    return state


def record_partial_stage_success(
    project_dir: Path,
    stage_name: str,
    *,
    message: str,
    outputs: list[str] | None = None,
) -> dict[str, Any]:
    """Close a successful smoke-test attempt while leaving the full stage pending."""
    state = load_run_state(project_dir)
    stage = next((row for row in state["stages"] if row["name"] == stage_name), None)
    if stage is None or stage["status"] != "running" or not stage["attempts"]:
        raise ValueError(f"stage is not running: {stage_name}")
    stage["attempts"][-1].update({
        "finished_at": _now(),
        "outcome": "smoke_test_completed",
        "message": message,
        "outputs": outputs or [],
    })
    stage["status"] = "pending"
    save_run_state(project_dir, state)
    return state


def interrupt_stage(project_dir: Path, stage_name: str, *, message: str) -> dict[str, Any]:
    state = load_run_state(project_dir)
    stage = next((row for row in state["stages"] if row["name"] == stage_name), None)
    if stage is None or stage["status"] != "running" or not stage["attempts"]:
        raise ValueError(f"stage is not running: {stage_name}")
    stage["status"] = "interrupted"
    stage["attempts"][-1].update({
        "finished_at": _now(), "outcome": "interrupted", "message": message,
    })
    save_run_state(project_dir, state)
    return state


def status_lines(state: dict[str, Any]) -> list[str]:
    lines = ["SEQUESTRA workflow status"]
    for stage in state["stages"]:
        attempts = len(stage.get("attempts", []))
        lines.append(f"[{stage['status'].upper():11}] {stage['name']} (attempts: {attempts})")
    next_stage = next_resumable_stage(state)
    lines.append(f"Next resumable stage: {next_stage or 'none; workflow complete'}")
    return lines
