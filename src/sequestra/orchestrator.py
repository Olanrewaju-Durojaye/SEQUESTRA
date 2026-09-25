"""Beginner-facing end-to-end workflow orchestration."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Callable

from .boltz2_affinity import run_boltz2_affinity
from .boltzgen_runner import run_boltzgen
from .catalytic_liability import run_catalytic_liability
from .reporting import run_reports
from .run_state import load_run_state, next_resumable_stage, status_lines
from .shortlist import create_shortlist, load_ranked_candidates, recommendation_for_total, selection_count
from .structural_confidence import run_structural_confidence

InputFunction = Callable[[str], str]
OutputFunction = Callable[[str], None]


def find_project(start: Path | None = None) -> Path:
    current = (start or Path.cwd()).expanduser().resolve()
    candidates = (current, *current.parents)
    for candidate in candidates:
        if (candidate / "configuration/project_specification.yaml").is_file():
            return candidate
    raise FileNotFoundError(
        "No SEQUESTRA project was found in this directory or its parents. "
        "Enter the project directory and run 'sequestra run' again."
    )


def find_config(explicit: Path | None = None) -> Path:
    candidates = []
    if explicit is not None:
        candidates.append(explicit.expanduser())
    environment = os.environ.get("SEQUESTRA_CONFIG")
    if environment:
        candidates.append(Path(environment).expanduser())
    candidates.extend((
        Path.home() / ".config/sequestra/sequestra.config.json",
        Path(__file__).resolve().parents[2] / "sequestra.config.json",
    ))
    for candidate in candidates:
        if candidate.is_file():
            resolved = candidate.resolve()
            if explicit is not None or environment:
                return resolved
            try:
                import json
                data = json.loads(resolved.read_text(encoding="utf-8"))
                conda = Path(data["conda_executable"]).expanduser()
                environments = data["environments"]
                if conda.is_file() and all(Path(value).expanduser().is_dir() for value in environments.values()):
                    return resolved
            except (KeyError, TypeError, ValueError, OSError):
                pass
    raise FileNotFoundError(
        "SEQUESTRA could not locate a valid configuration for this computer. "
        "Run 'sequestra setup' once."
    )


def _percentage(
    root: Path,
    supplied: float | None,
    *,
    ask: InputFunction,
    emit: OutputFunction,
) -> float:
    ranked, metric, _ = load_ranked_candidates(root)
    total = len(ranked)
    recommendation = recommendation_for_total(total)
    emit(f"\nBoltzGen ranking is ready: {total} candidates.")
    emit(f"Ranking metric: {metric} (higher is better).")
    emit(
        "Suggested percentage: "
        f"{recommendation['minimum_percent']}-{recommendation['maximum_percent']}%."
    )
    if supplied is not None:
        selection_count(total, supplied)
        return supplied
    while True:
        raw = ask("Percentage of top-ranked candidates to carry forward: ").strip()
        try:
            value = float(raw)
            count = selection_count(total, value)
        except ValueError as error:
            emit(str(error))
            continue
        confirmation = ask(f"Carry forward top {value:g}% ({count} candidates)? [Y/n]: ").strip().lower()
        if confirmation in {"", "y", "yes"}:
            return value
        emit("Choose the percentage again.")


def run_workflow(
    project_dir: Path | None = None,
    *,
    config_path: Path | None = None,
    devices: int = 1,
    shortlist_percentage: float | None = None,
    ask: InputFunction = input,
    emit: OutputFunction = print,
) -> int:
    if devices < 1:
        raise ValueError("devices must be at least 1")
    root = find_project(project_dir) if project_dir is not None else find_project()
    config = find_config(config_path)
    emit(f"\nSEQUESTRA project: {root}")
    emit("The workflow will resume from its first incomplete stage.")
    while True:
        state = load_run_state(root)
        stage = next_resumable_stage(state)
        if stage is None:
            emit("\n" + "\n".join(status_lines(state)))
            emit("Workflow complete.")
            return 0
        emit(f"\nStarting stage: {stage}")
        if stage == "boltzgen":
            code = run_boltzgen(root, config_path=config, devices=devices, emit=emit)
        elif stage == "shortlist":
            percentage = _percentage(root, shortlist_percentage, ask=ask, emit=emit)
            result = create_shortlist(root, percentage)
            emit(f"Shortlist completed: {result['selected']} of {result['total']} candidates retained.")
            code = 0
        elif stage == "catalytic_liability":
            code = run_catalytic_liability(root, sequestra_config_path=config)
        elif stage == "boltz2_affinity":
            code = run_boltz2_affinity(root, config_path=config)
        elif stage == "structural_confidence":
            code = run_structural_confidence(root, config_path=config)
        elif stage == "reports":
            code = run_reports(root)
        else:
            raise ValueError(f"unsupported workflow stage: {stage}")
        if code != 0:
            emit(f"Workflow paused at {stage}. Run 'sequestra run' again after resolving the reported issue.")
            return code
