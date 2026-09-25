"""Affinity-probability-only ranking and resumable shortlist materialization."""

from __future__ import annotations

import csv
import hashlib
import json
from math import ceil, isfinite
import os
from pathlib import Path
import shutil
import tempfile
from typing import Any
from uuid import uuid4

from .decision_rules import shortlist_recommendation
from .run_state import finish_stage, start_stage

AFFINITY_COLUMNS = (
    "affinity_probability_binary",
    "affinity_probability_binary1",
    "affinity_probability",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def metrics_path(project_dir: Path) -> Path:
    path = (
        project_dir.expanduser().resolve()
        / "boltzgen/run/final_ranked_designs/all_designs_metrics.csv"
    )
    if not path.is_file():
        raise FileNotFoundError(f"BoltzGen metrics file not found: {path}")
    return path


def load_ranked_candidates(project_dir: Path) -> tuple[list[dict[str, str]], str, Path]:
    source = metrics_path(project_dir)
    with source.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        metric = next((name for name in AFFINITY_COLUMNS if name in fieldnames), None)
        if metric is None:
            raise ValueError(
                "BoltzGen metrics do not contain a supported affinity-probability column; "
                f"observed columns: {', '.join(fieldnames)}"
            )
        rows = list(reader)
    if not rows:
        raise ValueError("BoltzGen metrics file contains no candidates")
    seen: set[str] = set()
    for index, row in enumerate(rows, start=1):
        candidate_id = (row.get("id") or "").strip()
        if not candidate_id:
            filename = (row.get("file_name") or "").strip()
            candidate_id = Path(filename).stem if filename else f"row_{index}"
        if candidate_id in seen:
            raise ValueError(f"duplicate BoltzGen candidate identifier: {candidate_id}")
        seen.add(candidate_id)
        row["sequestra_candidate_id"] = candidate_id
        try:
            value = float(row[metric])
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid {metric} value for {candidate_id}") from error
        if not isfinite(value):
            raise ValueError(f"non-finite {metric} value for {candidate_id}")
        row["sequestra_affinity_probability"] = repr(value)
    ranked = sorted(
        rows,
        key=lambda row: (
            -float(row["sequestra_affinity_probability"]),
            row["sequestra_candidate_id"],
        ),
    )
    for rank, row in enumerate(ranked, start=1):
        row["sequestra_rank"] = str(rank)
    return ranked, metric, source


def selection_count(total: int, percentage: float) -> int:
    if total < 1:
        raise ValueError("candidate total must be positive")
    if not isfinite(percentage) or not 0 < percentage <= 100:
        raise ValueError("percentage must be greater than 0 and at most 100")
    return max(1, ceil(total * percentage / 100.0))


def _write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _copy_candidate_structures(root: Path, selected: list[dict[str, str]], staging: Path) -> list[str]:
    inverse = root / "boltzgen/run/intermediate_designs_inverse_folded"
    copied: list[str] = []
    for row in selected:
        filename = (row.get("file_name") or "").strip()
        if not filename:
            raise ValueError(
                f"candidate {row['sequestra_candidate_id']} has no file_name for structure retrieval"
            )
        for source_subdir, destination_subdir, required in (
            ("refold_cif", "selected_complexes", True),
            ("refold_design_cif", "selected_binders", False),
        ):
            source = inverse / source_subdir / filename
            if not source.is_file():
                if required:
                    raise FileNotFoundError(f"selected candidate structure not found: {source}")
                continue
            destination = staging / destination_subdir / filename
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, destination)
            copied.append(str(destination.relative_to(staging)))
    return copied


def _update_project_specification(root: Path, percentage: float, count: int) -> None:
    path = root / "configuration/project_specification.yaml"
    specification = json.loads(path.read_text(encoding="utf-8"))
    rules = specification["decision_rules"]
    rules["top_fraction"] = percentage / 100.0
    rules["top_fraction_percent"] = percentage
    rules["top_fraction_selection_status"] = "selected_after_boltzgen_ranking"
    rules["selected_candidate_count"] = count
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(specification, handle, indent=2)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def _existing_selection_matches(destination: Path, percentage: float) -> bool:
    manifest_path = destination / "selection_manifest.json"
    if not manifest_path.is_file():
        return False
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if float(manifest["selected_percentage"]) != percentage:
        raise ValueError(
            "a validated shortlist already exists with a different percentage; "
            "create a new project to change a completed scientific decision"
        )
    required = (destination / "complete_affinity_ranking.csv", destination / "selected_candidates.csv")
    if not all(path.is_file() for path in required):
        raise ValueError("existing shortlist is incomplete")
    return True


def create_shortlist(project_dir: Path, percentage: float) -> dict[str, Any]:
    root = project_dir.expanduser().resolve()
    ranked, metric, source = load_ranked_candidates(root)
    count = selection_count(len(ranked), percentage)
    destination = root / "shortlist/selection"
    if destination.exists():
        if not _existing_selection_matches(destination, percentage):
            raise ValueError(f"refusing to overwrite unexpected shortlist directory: {destination}")
        start_stage(root, "shortlist", command=["internal:shortlist", str(percentage)])
        _update_project_specification(root, percentage, count)
        outputs = [
            "shortlist/selection/complete_affinity_ranking.csv",
            "shortlist/selection/selected_candidates.csv",
            "shortlist/selection/selection_manifest.json",
        ]
        finish_stage(root, "shortlist", succeeded=True, message="existing shortlist validated", outputs=outputs)
        return {"metric": metric, "total": len(ranked), "selected": count, "percentage": percentage}

    start_stage(root, "shortlist", command=["internal:shortlist", str(percentage)])
    staging = root / "shortlist" / f".selection-{uuid4().hex}.tmp"
    staging.mkdir(parents=True)
    try:
        selected = ranked[:count]
        full_fields = list(ranked[0].keys())
        if "sequestra_selected" not in full_fields:
            full_fields.append("sequestra_selected")
        for row in ranked:
            row["sequestra_selected"] = "true" if int(row["sequestra_rank"]) <= count else "false"
        _write_csv(staging / "complete_affinity_ranking.csv", ranked, full_fields)
        _write_csv(staging / "selected_candidates.csv", selected, full_fields)
        copied = _copy_candidate_structures(root, selected, staging)
        manifest = {
            "schema_version": 1,
            "source_metrics": str(source.relative_to(root)),
            "source_metrics_sha256": _sha256(source),
            "ranking_metric": metric,
            "ranking_direction": "descending_higher_is_better",
            "other_metrics_used_for_ranking": [],
            "tie_breaker": "candidate_id_ascending",
            "tie_boundary_policy": "retain_exact_ceiling_count_without_expanding_ties",
            "selected_percentage": percentage,
            "total_candidates": len(ranked),
            "selected_count": count,
            "selected_candidate_ids": [row["sequestra_candidate_id"] for row in selected],
            "copied_structures": copied,
        }
        (staging / "selection_manifest.json").write_text(
            json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
        )
        os.replace(staging, destination)
        _update_project_specification(root, percentage, count)
        outputs = [str(path.relative_to(root)) for path in destination.rglob("*") if path.is_file()]
        finish_stage(
            root, "shortlist", succeeded=True,
            message=f"selected top {percentage}% ({count} of {len(ranked)}) by {metric} only",
            outputs=outputs,
        )
        return {
            "metric": metric, "total": len(ranked), "selected": count,
            "percentage": percentage, "destination": str(destination),
        }
    except Exception as error:
        finish_stage(root, "shortlist", succeeded=False, message=str(error))
        raise


def recommendation_for_total(total: int) -> dict[str, Any]:
    return shortlist_recommendation(total)
