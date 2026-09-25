from __future__ import annotations

import csv
import json
import tempfile
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path

from sequestra.run_state import initialize_run_state, load_run_state, save_run_state
from sequestra.structural_confidence import run_structural_confidence


class StructuralConfidenceTests(unittest.TestCase):
    def project(self, base: Path, mode: str, ids: list[str]) -> tuple[Path, Path]:
        root = base / "project"
        (root / "configuration").mkdir(parents=True)
        (root / "boltz2_affinity/results").mkdir(parents=True)
        (root / "configuration/project_specification.yaml").write_text(
            json.dumps({"project_name": "test", "reference": {"mode": mode}})
        )
        (root / "boltz2_affinity/results/summary.json").write_text(
            json.dumps({"candidate_count": 2, "qualified_count": len(ids), "qualified_ids": ids})
        )
        (root / "boltz2_affinity/results/affinity_qualified_candidates.fasta").write_text(
            "".join(f">{identifier}\nACDEFG\n" for identifier in ids)
        )
        initialize_run_state(
            root,
            reference_mode=mode,
            catalytic_liability_enabled=mode == "catalytic-reference",
        )
        state = load_run_state(root)
        for stage in state["stages"]:
            if stage["name"] in {"boltzgen", "shortlist", "catalytic_liability", "boltz2_affinity"}:
                stage["status"] = "completed"
        save_run_state(root, state)
        config = base / "config.json"
        config.write_text(json.dumps({"structural_confidence": {
            "minimum_confidence_score": 0.7,
            "minimum_ligand_iptm": 0.6,
            "minimum_complex_iplddt": 0.7,
        }}))
        return root, config

    def artifacts(self, root: Path, identifier: str, confidence: float, ligand: float, interface: float) -> None:
        destination = root / f"boltz2_affinity/raw/{identifier}/predictions/{identifier}"
        destination.mkdir(parents=True)
        (destination / f"confidence_{identifier}_model_0.json").write_text(json.dumps({
            "confidence_score": confidence,
            "ptm": 0.8,
            "iptm": 0.75,
            "ligand_iptm": ligand,
            "protein_iptm": 0.0,
            "complex_plddt": 0.82,
            "complex_iplddt": interface,
            "complex_pde": 2.0,
            "complex_ipde": 3.0,
        }))
        (destination / f"{identifier}_model_0.cif").write_text("data_test\n")

    def test_reuses_outputs_and_applies_all_thresholds_in_both_modes(self) -> None:
        for mode in ("catalytic-reference", "non-catalytic-reference"):
            with self.subTest(mode=mode), tempfile.TemporaryDirectory() as directory:
                root, config = self.project(Path(directory), mode, ["design_a", "design_b"])
                self.artifacts(root, "design_a", 0.80, 0.70, 0.75)
                self.artifacts(root, "design_b", 0.90, 0.50, 0.85)
                self.assertEqual(run_structural_confidence(root, config_path=config), 0)
                table = root / "structural_confidence/results/structural_confidence.csv"
                with table.open() as handle:
                    rows = list(csv.DictReader(handle))
                self.assertEqual(rows[0]["decision"], "advance_to_reports")
                self.assertEqual(rows[1]["decision"], "does_not_advance")
                self.assertIn("ligand_iptm_below_0.6", rows[1]["detail"])
                self.assertEqual(
                    (root / "structural_confidence/results/confidence_qualified_candidates.fasta").read_text(),
                    ">design_a\nACDEFG\n",
                )
                plot = root / "structural_confidence/results/structural_confidence.svg"
                ET.parse(plot)
                svg = plot.read_text()
                self.assertNotIn(">design_a</text>", svg)
                self.assertIn("<title>design_a:", svg)
                state = load_run_state(root)
                self.assertEqual(next(s for s in state["stages"] if s["name"] == "structural_confidence")["status"], "completed")
                self.assertEqual(next(s for s in state["stages"] if s["name"] == "reports")["status"], "pending")

    def test_zero_survivors_skip_confidence_and_complete_terminal_report(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root, config = self.project(Path(directory), "non-catalytic-reference", [])
            self.assertEqual(run_structural_confidence(root, config_path=config), 0)
            state = load_run_state(root)
            structural = next(s for s in state["stages"] if s["name"] == "structural_confidence")
            reports = next(s for s in state["stages"] if s["name"] == "reports")
            self.assertEqual(structural["status"], "skipped")
            self.assertEqual(structural["skip_reason"], "no_affinity_qualified_candidates")
            self.assertEqual(reports["status"], "completed")
            self.assertIn(
                "computational negative result",
                (root / "reports/terminal_summary.md").read_text(),
            )
            self.assertTrue((root / "reports/reproducibility_manifest.json").is_file())

    def test_completed_candidate_checkpoint_is_reused(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root, config = self.project(Path(directory), "catalytic-reference", ["design_a"])
            self.artifacts(root, "design_a", 0.80, 0.70, 0.75)
            checkpoints = root / "structural_confidence/checkpoints"
            checkpoints.mkdir(parents=True)
            cached_row = {
                "candidate_id": "design_a", "confidence_score": 0.8, "ptm": 0.8,
                "iptm": 0.75, "ligand_iptm": 0.7, "protein_iptm": 0.0,
                "complex_plddt": 0.82, "complex_iplddt": 0.75, "complex_pde": 2.0,
                "complex_ipde": 3.0, "minimum_confidence_score": 0.7,
                "minimum_ligand_iptm": 0.6, "minimum_complex_iplddt": 0.7,
                "decision": "advance_to_reports", "detail": "all_required_confidence_thresholds_met",
                "source_confidence_json": "cached", "source_structure": "cached",
            }
            (checkpoints / "design_a.json").write_text(json.dumps({
                "schema_version": 1, "status": "completed",
                "thresholds": {"minimum_confidence_score": 0.7, "minimum_ligand_iptm": 0.6, "minimum_complex_iplddt": 0.7},
                "row": cached_row,
            }))
            confidence = next((root / "boltz2_affinity/raw/design_a").rglob("confidence_*.json"))
            confidence.unlink()
            self.assertEqual(run_structural_confidence(root, config_path=config), 0)
            self.assertIn("cached", (root / "structural_confidence/results/structural_confidence.csv").read_text())


if __name__ == "__main__":
    unittest.main()
