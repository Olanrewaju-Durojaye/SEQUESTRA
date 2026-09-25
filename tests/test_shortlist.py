from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from sequestra.run_state import finish_stage, initialize_run_state, load_run_state, start_stage
from sequestra.shortlist import create_shortlist, load_ranked_candidates, selection_count


class ShortlistTests(unittest.TestCase):
    def _project(self, root: Path) -> Path:
        project = root / "project"
        metrics_dir = project / "boltzgen/run/final_ranked_designs"
        complex_dir = project / "boltzgen/run/intermediate_designs_inverse_folded/refold_cif"
        binder_dir = project / "boltzgen/run/intermediate_designs_inverse_folded/refold_design_cif"
        for path in (metrics_dir, complex_dir, binder_dir, project / "configuration", project / "shortlist"):
            path.mkdir(parents=True, exist_ok=True)
        rows = [
            {"id": "design_c", "file_name": "c.cif", "affinity_probability_binary1": "0.2"},
            {"id": "design_b", "file_name": "b.cif", "affinity_probability_binary1": "0.9"},
            {"id": "design_a", "file_name": "a.cif", "affinity_probability_binary1": "0.9"},
        ]
        with (metrics_dir / "all_designs_metrics.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        for filename in ("a.cif", "b.cif", "c.cif"):
            (complex_dir / filename).write_text("complex")
            (binder_dir / filename).write_text("binder")
        specification = {
            "decision_rules": {
                "top_fraction": None,
                "top_fraction_selection_status": "pending_after_boltzgen_ranking",
            }
        }
        (project / "configuration/project_specification.yaml").write_text(
            json.dumps(specification)
        )
        initialize_run_state(
            project, reference_mode="catalytic-reference", catalytic_liability_enabled=True
        )
        start_stage(project, "boltzgen")
        finish_stage(project, "boltzgen", succeeded=True)
        return project

    def test_ranking_uses_affinity_only_and_candidate_id_tie_break(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project = self._project(Path(directory))
            ranked, metric, _ = load_ranked_candidates(project)
            self.assertEqual(metric, "affinity_probability_binary1")
            self.assertEqual(
                [row["sequestra_candidate_id"] for row in ranked],
                ["design_a", "design_b", "design_c"],
            )

    def test_percentage_uses_ceiling_rounding(self) -> None:
        self.assertEqual(selection_count(3, 10), 1)
        self.assertEqual(selection_count(11, 10), 2)

    def test_materializes_and_checkpoints_selection(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project = self._project(Path(directory))
            result = create_shortlist(project, 50)
            self.assertEqual(result["selected"], 2)
            selection = project / "shortlist/selection"
            self.assertTrue((selection / "selected_complexes/a.cif").is_file())
            self.assertTrue((selection / "selected_complexes/b.cif").is_file())
            self.assertTrue((selection / "selected_binders/a.cif").is_file())
            manifest = json.loads((selection / "selection_manifest.json").read_text())
            self.assertEqual(manifest["other_metrics_used_for_ranking"], [])
            self.assertEqual(manifest["selected_candidate_ids"], ["design_a", "design_b"])
            state = load_run_state(project)
            self.assertEqual(state["stages"][1]["status"], "completed")
            specification = json.loads(
                (project / "configuration/project_specification.yaml").read_text()
            )
            self.assertEqual(specification["decision_rules"]["top_fraction"], 0.5)

    def test_invalid_percentage_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            selection_count(10, 0)
        with self.assertRaises(ValueError):
            selection_count(10, 101)


if __name__ == "__main__":
    unittest.main()
