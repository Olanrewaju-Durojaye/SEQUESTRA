"""Guided, mode-aware SEQUESTRA project creation."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Callable

from .decision_rules import (
    CATALYTIC_REFERENCE, NON_CATALYTIC_REFERENCE, shortlist_recommendation,
)
from .project import initialize_reference_project

InputFunction = Callable[[str], str]
OutputFunction = Callable[[str], None]


def _required(prompt: str, ask: InputFunction) -> str:
    while True:
        value = ask(prompt).strip()
        if value:
            return value
        print("A value is required.")


def _new_project_directory(ask: InputFunction) -> str:
    """Reject an existing nonempty destination before collecting all other inputs."""
    while True:
        value = _required("New project directory: ", ask)
        destination = Path(value).expanduser()
        if destination.exists() and not destination.is_dir():
            print("That path exists and is not a directory. Choose a new project directory.")
            continue
        if destination.exists() and any(destination.iterdir()):
            print("That directory is not empty. Choose a new subdirectory for this project.")
            continue
        return value


def _integer(prompt: str, ask: InputFunction, *, minimum: int = 1) -> int:
    while True:
        raw = ask(prompt).strip()
        try:
            value = int(raw)
        except ValueError:
            print("Enter a whole number.")
            continue
        if value >= minimum:
            return value
        print(f"Enter a value of at least {minimum}.")


def _float_default(prompt: str, ask: InputFunction, default: float) -> float:
    while True:
        raw = ask(f"{prompt} [{default}]: ").strip()
        if not raw:
            return default
        try:
            value = float(raw)
        except ValueError:
            print("Enter a numeric value.")
            continue
        if value >= 0:
            return value
        print("Enter a nonnegative value.")


def _yes_no(prompt: str, ask: InputFunction, *, default: bool = False) -> bool:
    suffix = " [Y/n]: " if default else " [y/N]: "
    while True:
        raw = ask(prompt + suffix).strip().lower()
        if not raw:
            return default
        if raw in {"y", "yes"}:
            return True
        if raw in {"n", "no"}:
            return False
        print("Enter yes or no.")


def collect_answers(ask: InputFunction = input) -> dict[str, Any]:
    """Collect initialization values without collecting the post-ranking fraction."""
    print("\nSEQUESTRA guided project creation")
    print("1. Catalytic reference")
    print("2. Non-catalytic binding reference")
    while True:
        selection = ask("Reference mode [1/2]: ").strip()
        if selection in {"1", "2"}:
            break
        print("Choose 1 or 2.")
    mode = CATALYTIC_REFERENCE if selection == "1" else NON_CATALYTIC_REFERENCE
    answers: dict[str, Any] = {
        "project_dir": _new_project_directory(ask),
        "project_name": _required("Project name: ", ask),
        "reference_mode": mode,
        "reference_pdb": _required("Reference PDB file: ", ask),
        "reference_chain": _required("Reference PDB chain: ", ask),
        "complete_fasta": _required("Complete reference FASTA file: ", ask),
        "ligand_name": _required("Ligand name: ", ask),
        "ligand_ccd": _required("Ligand CCD: ", ask),
        "ligand_smiles": _required("Ligand isomeric SMILES: ", ask),
        "ligand_inchikey": ask("Ligand InChIKey (optional): ").strip() or None,
        "design_min_length": _integer("Minimum design length: ", ask, minimum=20),
        "design_max_length": _integer("Maximum design length: ", ask, minimum=20),
        "number_of_designs": _integer("Number of BoltzGen designs: ", ask),
        "affinity_comparison_margin": _float_default(
            "Boltz-2 reference comparison margin", ask, 0.0
        ),
        "top_fraction": None,
    }
    answers["catalytic_liability_enabled"] = (
        True if mode == CATALYTIC_REFERENCE else _yes_no(
            "Enable the optional prediction-only catalytic-liability screen?", ask
        )
    )
    if answers["catalytic_liability_enabled"]:
        answers["kcat_unit"] = ask("kcat unit [s^-1]: ").strip() or "s^-1"
        answers["km_unit"] = ask("Km unit [M]: ").strip() or "M"
    else:
        answers["kcat_unit"] = None
        answers["km_unit"] = None
    return answers


def load_answers(path: Path) -> dict[str, Any]:
    data = json.loads(path.expanduser().read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("answers file must contain one JSON object")
    return data


def _normalize_answers(answers: dict[str, Any]) -> dict[str, Any]:
    normalized = dict(answers)
    required = (
        "project_dir", "project_name", "reference_mode", "reference_pdb",
        "reference_chain", "complete_fasta", "ligand_name", "ligand_ccd",
        "ligand_smiles", "design_min_length", "design_max_length", "number_of_designs",
    )
    missing = [key for key in required if key not in normalized]
    if missing:
        raise ValueError("missing required answers: " + ", ".join(missing))
    for key in ("project_dir", "reference_pdb", "complete_fasta"):
        normalized[key] = Path(normalized[key])
    normalized.setdefault("ligand_inchikey", None)
    normalized.setdefault("top_fraction", None)
    normalized.setdefault("affinity_comparison_margin", 0.0)
    normalized.setdefault("catalytic_liability_enabled", False)
    normalized.setdefault("kcat_unit", None)
    normalized.setdefault("km_unit", None)
    return normalized


def print_summary(answers: dict[str, Any], emit: OutputFunction = print) -> None:
    emit("\nProject summary")
    fields = (
        ("Name", "project_name"), ("Mode", "reference_mode"),
        ("Destination", "project_dir"), ("Reference PDB", "reference_pdb"),
        ("Reference chain", "reference_chain"), ("Complete FASTA", "complete_fasta"),
        ("Ligand", "ligand_name"), ("CCD", "ligand_ccd"),
        ("Design length", None), ("Number of designs", "number_of_designs"),
    )
    for label, key in fields:
        value = (
            f"{answers['design_min_length']}-{answers['design_max_length']} residues"
            if key is None else answers[key]
        )
        emit(f"  {label}: {value}")
    recommendation = shortlist_recommendation(int(answers["number_of_designs"]))
    emit(
        "  Non-binding post-ranking recommendation: "
        f"{recommendation['minimum_percent']}-{recommendation['maximum_percent']}%"
    )
    emit("  Final shortlist percentage: pending user selection after BoltzGen ranking")
    liability_enabled = bool(answers.get("catalytic_liability_enabled"))
    emit(f"  Catalytic-liability screen: {'enabled' if liability_enabled else 'disabled'}")


def create_project(
    *, answers: dict[str, Any], assume_yes: bool = False,
    ask: InputFunction = input, emit: OutputFunction = print,
) -> dict[str, Any] | None:
    normalized = _normalize_answers(answers)
    print_summary(normalized, emit)
    if not assume_yes and not _yes_no("Create this project?", ask):
        emit("Project creation cancelled; no files were written.")
        return None
    result = initialize_reference_project(**normalized)
    recommendation = shortlist_recommendation(int(normalized["number_of_designs"]))
    result["shortlist_recommendation"] = recommendation
    emit("\nSEQUESTRA project created successfully.")
    emit(f"Project: {result['project_dir']}")
    emit("Shortlist selection remains pending until BoltzGen ranking is available.")
    emit(
        "Suggested starting range: "
        f"{recommendation['minimum_percent']}-{recommendation['maximum_percent']}%."
    )
    return result
