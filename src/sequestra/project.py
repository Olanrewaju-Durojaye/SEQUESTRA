"""Create validated, self-contained SEQUESTRA reference projects."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

from .decision_rules import (
    CATALYTIC_REFERENCE,
    NON_CATALYTIC_REFERENCE,
    DecisionSettings,
    mode_workflow,
    shortlist_recommendation,
)
from .reference import (
    expression_tag_residues, inspect_pdb, map_resolved_to_complete,
    read_fasta, select_complete_sequence,
)
from .run_state import initialize_run_state


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _validate_smiles(smiles: str) -> str:
    normalized = smiles.strip()
    if not normalized or any(char.isspace() for char in normalized):
        raise ValueError("ligand SMILES must be nonempty and contain no whitespace")
    pairs = {"(": ")", "[": "]"}
    stack: list[str] = []
    for char in normalized:
        if char in pairs:
            stack.append(pairs[char])
        elif char in pairs.values():
            if not stack or stack.pop() != char:
                raise ValueError("ligand SMILES contains unbalanced brackets")
    if stack:
        raise ValueError("ligand SMILES contains unbalanced brackets")
    return normalized


def _validate_design_settings(minimum: int, maximum: int, count: int) -> None:
    if minimum < 20:
        raise ValueError("minimum design length must be at least 20 residues")
    if maximum < minimum:
        raise ValueError("maximum design length cannot be less than minimum")
    if count < 1:
        raise ValueError("number of designs must be positive")


def initialize_reference_project(
    *,
    project_dir: Path,
    project_name: str,
    reference_mode: str,
    reference_pdb: Path,
    reference_chain: str,
    complete_fasta: Path,
    ligand_name: str,
    ligand_ccd: str,
    ligand_smiles: str,
    ligand_inchikey: str | None,
    design_min_length: int,
    design_max_length: int,
    number_of_designs: int,
    top_fraction: float | None = None,
    affinity_comparison_margin: float = 0.0,
    catalytic_liability_enabled: bool = False,
    kcat_unit: str | None = None,
    km_unit: str | None = None,
) -> dict[str, Any]:
    destination = project_dir.expanduser().resolve()
    if destination.exists() and any(destination.iterdir()):
        raise FileExistsError(f"refusing to overwrite nonempty project directory: {destination}")
    source_pdb = reference_pdb.expanduser().resolve()
    source_fasta = complete_fasta.expanduser().resolve()
    if not source_pdb.is_file() or not source_fasta.is_file():
        raise FileNotFoundError("reference PDB and complete FASTA must both exist")
    if len(reference_chain) != 1 or not reference_chain.strip():
        raise ValueError("reference chain must be one nonblank PDB chain identifier")
    settings = DecisionSettings(
        top_fraction=top_fraction,
        affinity_comparison_margin=affinity_comparison_margin,
        catalytic_liability_enabled=catalytic_liability_enabled,
    )
    settings.validate(reference_mode)
    kinetic_screen_enabled = reference_mode == CATALYTIC_REFERENCE or catalytic_liability_enabled
    if kinetic_screen_enabled:
        kcat_unit = (kcat_unit or "s^-1").strip()
        km_unit = (km_unit or "M").strip()
        if not kcat_unit or not km_unit:
            raise ValueError("kinetic units must be nonempty when the liability screen is enabled")
    else:
        kcat_unit = None
        km_unit = None
    _validate_design_settings(design_min_length, design_max_length, number_of_designs)
    ccd = ligand_ccd.strip().upper()
    if not ccd or len(ccd) > 5 or not ccd.isalnum():
        raise ValueError("ligand CCD must be a 1-5 character alphanumeric identifier")
    smiles = _validate_smiles(ligand_smiles)

    records = read_fasta(source_fasta)
    complete_record = select_complete_sequence(records, reference_chain)
    excluded_expression_tags = sorted(expression_tag_residues(source_pdb, reference_chain))
    residues, bound_ligands = inspect_pdb(source_pdb, reference_chain)
    if ccd not in bound_ligands:
        available = ", ".join(sorted(bound_ligands)) or "none"
        raise ValueError(
            f"ligand CCD {ccd} was not found as a bound non-water component on "
            f"chain {reference_chain}; observed: {available}"
        )
    mapping = map_resolved_to_complete(complete_record.sequence, residues)
    resolved_sequence = "".join(residue.amino_acid for residue in residues)

    inputs = destination / "inputs"
    configuration = destination / "configuration"
    stage_names = (
        "boltzgen", "shortlist", "catalytic_liability", "boltz2_affinity",
        "structural_confidence", "reports", "logs",
    )
    for directory in (inputs, configuration, *(destination / name for name in stage_names)):
        directory.mkdir(parents=True, exist_ok=True)
    pdb_copy = inputs / "reference_structure.pdb"
    fasta_copy = inputs / "reference_complete.fasta"
    resolved_fasta = inputs / "reference_observed.fasta"
    shutil.copy2(source_pdb, pdb_copy)
    shutil.copy2(source_fasta, fasta_copy)
    resolved_fasta.write_text(
        f">{project_name}|chain_{reference_chain}|coordinate_resolved\n"
        + "\n".join(resolved_sequence[i:i + 70] for i in range(0, len(resolved_sequence), 70))
        + "\n",
        encoding="utf-8",
    )

    mapping_path = inputs / "reference_residue_mapping.csv"
    with mapping_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping[0]))
        writer.writeheader()
        writer.writerows(mapping)

    ligand_record = {
        "name": ligand_name.strip(),
        "ccd": ccd,
        "smiles": smiles,
        "inchikey": (ligand_inchikey or "").strip(),
        "bound_component_verified": True,
        "pdb_chain": reference_chain,
    }
    ligand_path = inputs / "ligand_record.yaml"
    ligand_path.write_text(json.dumps(ligand_record, indent=2) + "\n", encoding="utf-8")
    modifications = {
        "schema_version": 1,
        "source_structure_preserved_unchanged": True,
        "coordinate_validation_exclusions": [
            {
                "chain": reference_chain,
                "pdb_author_number": number,
                "pdb_insertion_code": insertion,
                "reason": "PDB_SEQADV_EXPRESSION_TAG",
            }
            for number, insertion in excluded_expression_tags
        ],
        "note": (
            "Residues explicitly annotated as EXPRESSION TAG in PDB SEQADV records "
            "are excluded from canonical FASTA mapping; the copied PDB is preserved unchanged."
        ),
    }
    modifications_path = inputs / "reference_input_modifications.json"
    modifications_path.write_text(json.dumps(modifications, indent=2) + "\n", encoding="utf-8")
    unresolved_count = sum(not bool(row["resolved"]) for row in mapping)
    specification = {
        "schema_version": 2,
        "project_name": project_name.strip(),
        "reference": {
            "mode": reference_mode,
            "chain": reference_chain,
            "complete_fasta_header": complete_record.header,
            "complete_length": len(complete_record.sequence),
            "resolved_length": len(resolved_sequence),
            "unresolved_count": unresolved_count,
            "structure": str(pdb_copy.relative_to(destination)),
            "complete_fasta": str(fasta_copy.relative_to(destination)),
            "resolved_fasta": str(resolved_fasta.relative_to(destination)),
            "residue_mapping": str(mapping_path.relative_to(destination)),
            "input_modifications": str(modifications_path.relative_to(destination)),
            "excluded_expression_tag_residue_count": len(excluded_expression_tags),
        },
        "ligand": ligand_record,
        "design": {
            "minimum_length": design_min_length,
            "maximum_length": design_max_length,
            "number_of_designs": number_of_designs,
            "initial_ranking_metric_family": "boltzgen_affinity_probability",
            "accepted_ranking_columns": [
                "affinity_probability_binary",
                "affinity_probability_binary1",
                "affinity_probability",
            ],
            "initial_ranking_direction": "higher",
        },
        "decision_rules": {
            "top_fraction": top_fraction,
            "top_fraction_selection_status": "selected" if top_fraction is not None else "pending_after_boltzgen_ranking",
            "top_fraction_rounding": "ceiling",
            "shortlist_recommendation": shortlist_recommendation(number_of_designs),
            "affinity_comparison_margin": affinity_comparison_margin,
            "boltz2_affinity_direction": "higher_is_better",
            "kinetic_threshold_source": "reference_predictions",
            "kcat_adverse_condition": "candidate_kcat > reference_kcat",
            "km_adverse_condition": "candidate_km < reference_km",
            "exclusion_logic": "both_conditions_required",
            "kcat_unit": kcat_unit,
            "km_unit": km_unit,
            "catalytic_liability_enabled": (
                True if reference_mode == CATALYTIC_REFERENCE else catalytic_liability_enabled
            ),
            "workflow": mode_workflow(reference_mode, catalytic_liability_enabled),
        },
        "input_sha256": {
            "reference_structure": _sha256(pdb_copy),
            "reference_complete_fasta": _sha256(fasta_copy),
            "reference_resolved_fasta": _sha256(resolved_fasta),
            "reference_input_modifications": _sha256(modifications_path),
        },
        "scientific_interpretation": {
            "kinetic_predictions_establish_non_catalysis": False,
            "predicted_affinity_establishes_experimental_binding": False,
            "catalytic_liability_screen_is_prediction_only": True,
            "noncatalytic_mode_liability_result_is_not_a_catalytic_classification": True,
            "prediction_failure_establishes_non_catalysis": False,
        },
    }
    specification_path = configuration / "project_specification.yaml"
    specification_path.write_text(json.dumps(specification, indent=2) + "\n", encoding="utf-8")
    commands = {
        "schema_version": 1,
        "execution_policy": "sequential",
        "commands": {
            name: {"command": None, "status": "not_configured"}
            for name in (
                "boltzgen", "shortlist", "catalytic_liability",
                "boltz2_affinity", "structural_confidence",
            )
        },
    }
    (configuration / "tool_commands.json").write_text(
        json.dumps(commands, indent=2) + "\n", encoding="utf-8"
    )
    run_state = initialize_run_state(
        destination,
        reference_mode=reference_mode,
        catalytic_liability_enabled=catalytic_liability_enabled,
    )
    return {
        "project_dir": str(destination),
        "project_specification": str(specification_path),
        "complete_length": len(complete_record.sequence),
        "resolved_length": len(resolved_sequence),
        "unresolved_count": unresolved_count,
        "excluded_expression_tag_residue_count": len(excluded_expression_tags),
        "ligand_ccd": ccd,
        "reference_mode": reference_mode,
        "next_resumable_stage": next(
            row["name"] for row in run_state["stages"] if row["status"] == "pending"
        ),
    }
