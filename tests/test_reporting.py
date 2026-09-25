from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from sequestra.reporting import run_reports
from sequestra.run_state import initialize_run_state, load_run_state, save_run_state


class ReportingTests(unittest.TestCase):
    def test_survivor_report_completes_workflow_and_writes_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "project"
            (root / "configuration").mkdir(parents=True)
            (root / "configuration/project_specification.yaml").write_text(json.dumps({
                "project_name": "pilot", "reference": {"mode": "non-catalytic-reference"}
            }))
            initialize_run_state(root, reference_mode="non-catalytic-reference", catalytic_liability_enabled=False)
            state = load_run_state(root)
            for stage in state["stages"]:
                if stage["name"] in {"boltzgen", "shortlist", "boltz2_affinity", "structural_confidence"}:
                    stage["status"] = "completed"
            save_run_state(root, state)
            affinity = root / "boltz2_affinity/results"
            confidence = root / "structural_confidence/results"
            affinity.mkdir(parents=True)
            confidence.mkdir(parents=True)
            (affinity / "summary.json").write_text(json.dumps({"candidate_count": 10, "qualified_count": 1}))
            with (affinity / "boltz2_affinity.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["candidate_id", "affinity_probability_binary"])
                writer.writeheader(); writer.writerow({"candidate_id": "design_a", "affinity_probability_binary": 0.8})
            (confidence / "summary.json").write_text(json.dumps({"candidate_count": 1, "qualified_count": 1}))
            with (confidence / "structural_confidence.csv").open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=["candidate_id", "confidence_score", "decision"])
                writer.writeheader(); writer.writerow({"candidate_id": "design_a", "confidence_score": 0.9, "decision": "advance_to_reports"})
            (confidence / "confidence_qualified_candidates.fasta").write_text(">design_a\nACDE\n")
            self.assertEqual(run_reports(root), 0)
            self.assertIn("design_a", (root / "reports/final_candidates.csv").read_text())
            self.assertEqual((root / "reports/final_candidates.fasta").read_text(), ">design_a\nACDE\n")
            manifest = json.loads((root / "reports/reproducibility_manifest.json").read_text())
            self.assertTrue(any(row["path"] == "reports/workflow_summary.json" for row in manifest["files"]))
            self.assertFalse(any(row["path"] == "configuration/run_state.json" for row in manifest["files"]))
            self.assertEqual(next(s for s in load_run_state(root)["stages"] if s["name"] == "reports")["status"], "completed")


if __name__ == "__main__":
    unittest.main()
