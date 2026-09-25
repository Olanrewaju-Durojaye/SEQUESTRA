"""One-time guided workstation configuration."""
from __future__ import annotations

import json
from pathlib import Path
import shutil
import subprocess
from typing import Callable

InputFunction = Callable[[str], str]
OutputFunction = Callable[[str], None]


def _directory(label: str, candidates: list[Path], ask: InputFunction) -> Path:
    for candidate in candidates:
        if candidate.expanduser().is_dir():
            return candidate.expanduser().resolve()
    while True:
        value = ask(f"{label} directory: ").strip()
        path = Path(value).expanduser()
        if path.is_dir():
            return path.resolve()
        print("That directory was not found. Please check the path.")


def _conda_environments(conda: Path) -> dict[str, Path]:
    completed = subprocess.run(
        [str(conda), "env", "list", "--json"],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    if completed.returncode != 0:
        return {}
    data = json.loads(completed.stdout)
    return {Path(value).name: Path(value).resolve() for value in data.get("envs", [])}


def configure_workstation(
    *, ask: InputFunction = input, emit: OutputFunction = print
) -> Path:
    emit("\nSEQUESTRA one-time workstation setup")
    found_conda = shutil.which("conda")
    if found_conda:
        conda = Path(found_conda).resolve()
        emit(f"Conda detected: {conda}")
    else:
        while True:
            value = ask("Conda executable: ").strip()
            conda = Path(value).expanduser()
            if conda.is_file():
                conda = conda.resolve()
                break
            emit("That executable was not found.")
    discovered = _conda_environments(conda)
    environments = {}
    for name in ("sequestra", "boltzgen", "dlkcat", "catpred", "boltz2"):
        if name in discovered:
            environments[name] = str(discovered[name])
            emit(f"Environment detected: {name}")
        else:
            environments[name] = str(_directory(f"Conda environment '{name}'", [], ask))
    home = Path.home()
    dlkcat = _directory("DLKcat", [home / "Applications/DLKcat"], ask)
    catpred = _directory(
        "CatPred",
        [home / "Applications/catpred_pipeline/CatPred", home / "Applications/CatPred"],
        ask,
    )
    catpred_data = _directory(
        "CatPred data",
        [catpred.parent / "data", home / "Applications/catpred_pipeline/data"],
        ask,
    )
    package_root = Path(__file__).resolve().parents[2]
    payload = {
        "schema_version": 1,
        "conda_executable": str(conda),
        "environments": environments,
        "repositories": {"sequestra": str(package_root)},
        "predictors": {
            "dlkcat": {"repository": str(dlkcat), "adapter": "native"},
            "catpred": {
                "repository": str(catpred),
                "data_root": str(catpred_data),
                "adapter": "native",
                "use_gpu": True,
            },
        },
        "boltz2": {
            "devices": 1,
            "accelerator": "gpu",
            "use_msa_server": True,
            "use_potentials": False,
        },
        "boltzgen": {
            "num_workers": 0,
            "minimum_available_ram_gib": 6.0,
            "memory_poll_seconds": 2.0,
        },
        "structural_confidence": {
            "minimum_confidence_score": 0.7,
            "minimum_ligand_iptm": 0.6,
            "minimum_complex_iplddt": 0.7,
        },
        "minimum_resources": {
            "free_disk_gib": 40,
            "available_ram_gib": 12,
            "gpu_memory_mib": 16000,
        },
    }
    destination = home / ".config/sequestra/sequestra.config.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    emit(f"Configuration saved: {destination}")
    return destination
