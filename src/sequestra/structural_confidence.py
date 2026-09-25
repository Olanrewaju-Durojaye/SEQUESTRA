"""Structural-confidence gate reusing preserved Boltz-2 predictions."""
from __future__ import annotations

import csv
import html
import json
import math
from pathlib import Path
from typing import Any

from .reference import read_fasta
from .run_state import finish_stage, load_run_state, skip_stage, start_stage

METRICS = ("confidence_score", "ligand_iptm", "complex_iplddt")
DEFAULT_THRESHOLDS = {
    "minimum_confidence_score": 0.70,
    "minimum_ligand_iptm": 0.60,
    "minimum_complex_iplddt": 0.70,
}


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def _number(data: dict[str, Any], key: str, *, bounded: bool = False) -> float | None:
    value = data.get(key)
    if value is None:
        return None
    if not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"invalid {key}: {value!r}")
    result = float(value)
    if bounded and not 0 <= result <= 1:
        raise ValueError(f"{key} must be between 0 and 1: {result}")
    return result


def load_thresholds(config_path: Path) -> dict[str, float]:
    config = _load(config_path.expanduser().resolve())
    configured = config.get("structural_confidence", {})
    thresholds = {
        key: float(configured.get(key, default))
        for key, default in DEFAULT_THRESHOLDS.items()
    }
    for key, value in thresholds.items():
        if not math.isfinite(value) or not 0 <= value <= 1:
            raise ValueError(f"{key} must be a finite value between 0 and 1")
    return thresholds


def qualified_ids(project_dir: Path) -> list[str]:
    root = project_dir.expanduser().resolve()
    summary_path = root / "boltz2_affinity/results/summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(f"Boltz-2 affinity summary not found: {summary_path}")
    ids = _load(summary_path).get("qualified_ids")
    if not isinstance(ids, list) or any(not isinstance(item, str) or not item for item in ids):
        raise ValueError("Boltz-2 affinity summary has invalid qualified_ids")
    return ids


def _artifacts(root: Path, identifier: str) -> tuple[Path, Path]:
    base = root / f"boltz2_affinity/raw/{identifier}"
    confidence = list(base.rglob(f"confidence_{identifier}_model_0.json"))
    if len(confidence) != 1:
        raise FileNotFoundError(
            f"expected exactly one preserved Boltz-2 confidence JSON for {identifier}; found {len(confidence)}"
        )
    structures = list(base.rglob(f"{identifier}_model_0.cif"))
    structures += list(base.rglob(f"{identifier}_model_0.pdb"))
    if len(structures) != 1:
        raise FileNotFoundError(
            f"expected exactly one preserved Boltz-2 model_0 structure for {identifier}; found {len(structures)}"
        )
    return confidence[0], structures[0]


def _checkpoint(root: Path, identifier: str) -> Path:
    return root / f"structural_confidence/checkpoints/{identifier}.json"


def _evaluate(root: Path, identifier: str, thresholds: dict[str, float]) -> dict[str, Any]:
    confidence_path, structure_path = _artifacts(root, identifier)
    data = _load(confidence_path)
    values = {metric: _number(data, metric, bounded=True) for metric in METRICS}
    missing = [metric for metric, value in values.items() if value is None]
    if missing:
        raise ValueError(
            f"required Boltz-2 confidence metric(s) missing for {identifier}: {', '.join(missing)}"
        )
    failures = []
    for metric in METRICS:
        threshold_key = f"minimum_{metric}"
        if values[metric] < thresholds[threshold_key]:
            failures.append(f"{metric}_below_{thresholds[threshold_key]:g}")
    passed = not failures
    return {
        "candidate_id": identifier,
        "confidence_score": values["confidence_score"],
        "ptm": _number(data, "ptm", bounded=True),
        "iptm": _number(data, "iptm", bounded=True),
        "ligand_iptm": values["ligand_iptm"],
        "protein_iptm": _number(data, "protein_iptm", bounded=True),
        "complex_plddt": _number(data, "complex_plddt", bounded=True),
        "complex_iplddt": values["complex_iplddt"],
        "complex_pde": _number(data, "complex_pde"),
        "complex_ipde": _number(data, "complex_ipde"),
        **thresholds,
        "decision": "advance_to_reports" if passed else "does_not_advance",
        "detail": "all_required_confidence_thresholds_met" if passed else ";".join(failures),
        "source_confidence_json": str(confidence_path.relative_to(root)),
        "source_structure": str(structure_path.relative_to(root)),
    }


def _write_plot(path: Path, rows: list[dict[str, Any]], thresholds: dict[str, float]) -> None:
    left, right = 90, 850
    panel_top = (55, 260, 465)
    panel_height = 150
    bottom = panel_top[-1] + panel_height
    count = len(rows)
    xmax = max(1, count - 1)
    sx = lambda index: left + (0.5 if count == 1 else index / xmax) * (right - left)
    sy = lambda value, top: top + panel_height - value * panel_height
    labels = {
        "confidence_score": "Overall confidence score",
        "ligand_iptm": "Protein-ligand ipTM",
        "complex_iplddt": "Interface-weighted pLDDT",
    }
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="720" viewBox="0 0 900 720">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:sans-serif;font-size:12px}.title{font-size:18px;font-weight:bold}.tick{font-size:11px;fill:#444}.grid{stroke:#e0e0e0;stroke-width:1}.panel{font-size:13px;font-weight:bold}</style>',
        '<text class="title" x="90" y="29">SEQUESTRA Boltz-2 structural-confidence assessment</text>',
    ]
    for metric, top in zip(METRICS, panel_top):
        threshold = thresholds[f"minimum_{metric}"]
        parts.append(f'<rect x="{left}" y="{top}" width="{right-left}" height="{max(0, sy(threshold, top)-top):.1f}" fill="#e8f5e9"/>')
        for index in range(5):
            value = index / 4
            y = sy(value, top)
            parts += [
                f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}"/>',
                f'<text class="tick" x="{left-12}" y="{y+4:.1f}" text-anchor="end">{value:.2f}</text>',
            ]
        y_threshold = sy(threshold, top)
        parts += [
            f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+panel_height}" stroke="black"/>',
            f'<line x1="{left}" y1="{top+panel_height}" x2="{right}" y2="{top+panel_height}" stroke="black"/>',
            f'<line x1="{left}" y1="{y_threshold:.1f}" x2="{right}" y2="{y_threshold:.1f}" stroke="#2e7d32" stroke-width="2" stroke-dasharray="7 5"/>',
            f'<text class="panel" x="{left}" y="{top-8}">{labels[metric]}</text>',
            f'<text x="{right-4}" y="{max(top+14,y_threshold-6):.1f}" text-anchor="end" fill="#1b5e20">pass if ≥ {threshold:.2f}</text>',
        ]
        for index, row in enumerate(rows):
            color = "#2e7d32" if row["decision"] == "advance_to_reports" else "#757575"
            title = html.escape(
                f"{row['candidate_id']}: {metric}={float(row[metric]):.6f}; {row['decision']}"
            )
            parts.append(
                f'<circle cx="{sx(index):.1f}" cy="{sy(float(row[metric]),top):.1f}" r="6" fill="{color}"><title>{title}</title></circle>'
            )
    tick_step = max(1, math.ceil(count / 10))
    for index in range(count):
        if index == 0 or index == count - 1 or (index + 1) % tick_step == 0:
            x = sx(index)
            parts += [
                f'<line x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom+6}" stroke="black"/>',
                f'<text class="tick" x="{x:.1f}" y="{bottom+22}" text-anchor="middle">{index+1}</text>',
            ]
    parts += [
        f'<text x="{(left+right)/2:.1f}" y="{bottom+50}" text-anchor="middle">Affinity-qualified candidates (input order)</text>',
        '<circle cx="250" cy="698" r="5" fill="#2e7d32"/><text x="261" y="702">Advance</text>',
        '<circle cx="350" cy="698" r="5" fill="#757575"/><text x="361" y="702">Does not advance</text>',
        '</svg>',
    ]
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def _terminal_report(root: Path, affinity_summary: dict[str, Any]) -> dict[str, Path]:
    from .reporting import _manifest

    specification = _load(root / "configuration/project_specification.yaml")
    destination = root / "reports"
    destination.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema_version": 1,
        "project_name": specification.get("project_name"),
        "reference_mode": specification["reference"]["mode"],
        "terminal_stage": "boltz2_affinity",
        "outcome": "no_candidate_exceeded_reference_affinity_benchmark",
        "candidate_count": affinity_summary.get("candidate_count", 0),
        "affinity_qualified_count": 0,
        "structural_confidence_status": "skipped",
        "scientific_interpretation": "No evaluated design had affinity_probability_binary strictly greater than the reference plus the configured margin. This is a computational negative result and does not establish absence of experimental binding.",
    }
    json_path = destination / "terminal_summary.json"
    _write_json(json_path, payload)
    markdown_path = destination / "terminal_summary.md"
    markdown_path.write_text(
        "# SEQUESTRA terminal summary\n\n"
        f"- Project: {payload['project_name']}\n"
        f"- Reference mode: {payload['reference_mode']}\n"
        "- Terminal scientific stage: Boltz-2 affinity benchmarking\n"
        f"- Candidates evaluated: {payload['candidate_count']}\n"
        "- Candidates exceeding the reference-plus-margin benchmark: 0\n"
        "- Structural-confidence assessment: skipped (no eligible candidates)\n\n"
        "No evaluated design had `affinity_probability_binary` strictly greater than "
        "the reference plus the configured margin. This is a computational negative "
        "result and does not establish absence of experimental binding.\n",
        encoding="utf-8",
    )
    manifest_path = destination / "reproducibility_manifest.json"
    manifest_path.write_text(
        json.dumps({
            "schema_version": 1,
            "project_root_name": root.name,
            "files": _manifest(
                root,
                {manifest_path, root / "configuration/run_state.json"},
            ),
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    return {"json": json_path, "markdown": markdown_path, "manifest": manifest_path}


def finalize_zero_survivor_workflow(project_dir: Path) -> dict[str, Path]:
    root = project_dir.expanduser().resolve()
    summary = _load(root / "boltz2_affinity/results/summary.json")
    if int(summary.get("qualified_count", -1)) != 0:
        raise ValueError("zero-survivor finalization requires qualified_count = 0")
    state = load_run_state(root)
    structural = next(row for row in state["stages"] if row["name"] == "structural_confidence")
    if structural["status"] not in {"skipped", "completed"}:
        skip_stage(
            root,
            "structural_confidence",
            reason="no_affinity_qualified_candidates",
            outputs=["boltz2_affinity/results/affinity_qualified_candidates.fasta"],
        )
    state = load_run_state(root)
    reports = next(row for row in state["stages"] if row["name"] == "reports")
    if reports["status"] == "completed":
        return {
            "json": root / "reports/terminal_summary.json",
            "markdown": root / "reports/terminal_summary.md",
            "manifest": root / "reports/reproducibility_manifest.json",
        }
    start_stage(root, "reports", command=["generate-terminal-zero-survivor-report"])
    try:
        paths = _terminal_report(root, summary)
        finish_stage(
            root,
            "reports",
            succeeded=True,
            message="terminal report generated: no candidates passed affinity benchmarking",
            outputs=[str(path.relative_to(root)) for path in paths.values()],
        )
        return paths
    except Exception as error:
        finish_stage(root, "reports", succeeded=False, message=str(error))
        raise


def run_structural_confidence(project_dir: Path, *, config_path: Path) -> int:
    root = project_dir.expanduser().resolve()
    ids = qualified_ids(root)
    if not ids:
        paths = finalize_zero_survivor_workflow(root)
        print("Structural confidence skipped: no candidates passed affinity benchmarking.")
        print(f"Terminal report: {paths['markdown']}")
        return 0
    thresholds = load_thresholds(config_path)
    start_stage(root, "structural_confidence", command=["reuse-preserved-boltz2-confidence", *ids])
    try:
        rows = []
        for identifier in ids:
            checkpoint = _checkpoint(root, identifier)
            if checkpoint.is_file():
                saved = _load(checkpoint)
                if (
                    saved.get("status") == "completed"
                    and saved.get("thresholds") == thresholds
                    and isinstance(saved.get("row"), dict)
                ):
                    rows.append(saved["row"])
                    continue
            row = _evaluate(root, identifier, thresholds)
            _write_json(
                checkpoint,
                {"schema_version": 1, "status": "completed", "thresholds": thresholds, "row": row},
            )
            rows.append(row)
        destination = root / "structural_confidence/results"
        destination.mkdir(parents=True, exist_ok=True)
        table = destination / "structural_confidence.csv"
        with table.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        passed = [row for row in rows if row["decision"] == "advance_to_reports"]
        sequences = {
            record.header.split()[0]: record.sequence
            for record in read_fasta(root / "boltz2_affinity/results/affinity_qualified_candidates.fasta")
        }
        fasta = destination / "confidence_qualified_candidates.fasta"
        fasta.write_text(
            "".join(f">{row['candidate_id']}\n{sequences[row['candidate_id']]}\n" for row in passed),
            encoding="utf-8",
        )
        plot = destination / "structural_confidence.svg"
        _write_plot(plot, rows, thresholds)
        summary = {
            "schema_version": 1,
            "source": "preserved_boltz2_model_0_outputs",
            "thresholds": thresholds,
            "candidate_count": len(rows),
            "qualified_count": len(passed),
            "qualified_ids": [row["candidate_id"] for row in passed],
            "interpretation_warning": "Confidence thresholds are configurable SEQUESTRA policy defaults, not official Boltz cutoffs or experimental validation.",
        }
        summary_path = destination / "summary.json"
        _write_json(summary_path, summary)
        outputs = [table, fasta, plot, summary_path]
        finish_stage(
            root,
            "structural_confidence",
            succeeded=True,
            message=f"{len(passed)} of {len(rows)} affinity-qualified candidates passed structural-confidence thresholds",
            outputs=[str(path.relative_to(root)) for path in outputs],
        )
        print(f"Structural confidence completed: {len(passed)} of {len(rows)} candidates advance to reports.")
        return 0
    except Exception as error:
        finish_stage(root, "structural_confidence", succeeded=False, message=str(error))
        raise
