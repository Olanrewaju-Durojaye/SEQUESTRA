"""Final human-readable and machine-readable workflow reporting."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from .run_state import finish_stage, load_run_state, start_stage


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _manifest(root: Path, excluded: set[Path]) -> list[dict[str, Any]]:
    records = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        if path in excluded or any(part.startswith(".") for part in path.relative_to(root).parts):
            continue
        records.append({
            "path": str(path.relative_to(root)),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        })
    return records


def run_reports(project_dir: Path) -> int:
    root = project_dir.expanduser().resolve()
    start_stage(root, "reports", command=["generate-final-workflow-report"])
    try:
        specification = _load(root / "configuration/project_specification.yaml")
        affinity = _load(root / "boltz2_affinity/results/summary.json")
        confidence = _load(root / "structural_confidence/results/summary.json")
        affinity_rows = {row["candidate_id"]: row for row in _rows(root / "boltz2_affinity/results/boltz2_affinity.csv")}
        confidence_rows = _rows(root / "structural_confidence/results/structural_confidence.csv")
        liability_rows = {row["candidate_id"]: row for row in _rows(root / "catalytic_liability/results/catalytic_liability.csv")}
        final = []
        for row in confidence_rows:
            if row["decision"] != "advance_to_reports":
                continue
            identifier = row["candidate_id"]
            combined = {"candidate_id": identifier}
            if identifier in affinity_rows:
                combined.update({f"affinity_{key}": value for key, value in affinity_rows[identifier].items() if key != "candidate_id"})
            combined.update({f"confidence_{key}": value for key, value in row.items() if key != "candidate_id"})
            if identifier in liability_rows:
                combined.update({f"liability_{key}": value for key, value in liability_rows[identifier].items() if key != "candidate_id"})
            final.append(combined)
        destination = root / "reports"
        destination.mkdir(parents=True, exist_ok=True)
        table = destination / "final_candidates.csv"
        if final:
            with table.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=list(final[0]))
                writer.writeheader()
                writer.writerows(final)
        else:
            table.write_text("candidate_id\n", encoding="utf-8")
        source_fasta = root / "structural_confidence/results/confidence_qualified_candidates.fasta"
        fasta = destination / "final_candidates.fasta"
        fasta.write_text(source_fasta.read_text(encoding="utf-8"), encoding="utf-8")
        outcome = "completed_with_final_candidates" if final else "no_candidate_passed_structural_confidence"
        summary = {
            "schema_version": 1,
            "project_name": specification.get("project_name"),
            "reference_mode": specification["reference"]["mode"],
            "outcome": outcome,
            "affinity_candidate_count": affinity["candidate_count"],
            "affinity_qualified_count": affinity["qualified_count"],
            "structural_confidence_candidate_count": confidence["candidate_count"],
            "final_candidate_count": len(final),
            "final_candidate_ids": [row["candidate_id"] for row in final],
            "scientific_interpretation": "All results are computational predictions requiring experimental validation.",
        }
        summary_json = destination / "workflow_summary.json"
        summary_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        summary_md = destination / "workflow_summary.md"
        summary_md.write_text(
            "# SEQUESTRA workflow summary\n\n"
            f"- Project: {summary['project_name']}\n"
            f"- Reference mode: {summary['reference_mode']}\n"
            f"- Outcome: {outcome}\n"
            f"- Candidates entering affinity benchmarking: {summary['affinity_candidate_count']}\n"
            f"- Candidates entering structural-confidence assessment: {summary['affinity_qualified_count']}\n"
            f"- Final computational candidates: {summary['final_candidate_count']}\n\n"
            "These results are computational predictions and require experimental validation.\n",
            encoding="utf-8",
        )
        manifest_path = destination / "reproducibility_manifest.json"
        manifest = {
            "schema_version": 1,
            "project_root_name": root.name,
            "files": _manifest(root, {manifest_path, root / "configuration/run_state.json"}),
        }
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
        outputs = [table, fasta, summary_json, summary_md, manifest_path]
        finish_stage(
            root,
            "reports",
            succeeded=True,
            message=f"final report generated with {len(final)} candidates",
            outputs=[str(path.relative_to(root)) for path in outputs],
        )
        print(f"Reports completed: {len(final)} final computational candidates.")
        print(f"Summary: {summary_md}")
        return 0
    except Exception as error:
        finish_stage(root, "reports", succeeded=False, message=str(error))
        raise
