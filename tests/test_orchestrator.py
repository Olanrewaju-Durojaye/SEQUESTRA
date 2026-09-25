from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from sequestra.orchestrator import find_project, run_workflow
from sequestra.run_state import finish_stage, initialize_run_state, load_run_state, start_stage


class OrchestratorTests(unittest.TestCase):
    def test_project_is_discovered_from_a_child_directory(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory) / "project"
            child = root / "boltzgen/run"
            child.mkdir(parents=True)
            (root / "configuration").mkdir()
            (root / "configuration/project_specification.yaml").write_text("{}")
            self.assertEqual(find_project(child), root.resolve())

    def test_single_run_command_executes_mode_aware_stage_sequence(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            root = base / "project"
            (root / "configuration").mkdir(parents=True)
            (root / "configuration/project_specification.yaml").write_text("{}")
            initialize_run_state(
                root,
                reference_mode="non-catalytic-reference",
                catalytic_liability_enabled=False,
            )
            config = base / "config.json"
            config.write_text("{}")
            calls: list[str] = []

            def complete(name: str):
                def runner(*args, **kwargs):
                    calls.append(name)
                    start_stage(root, name)
                    finish_stage(root, name, succeeded=True)
                    return 0
                return runner

            def shortlist(*args, **kwargs):
                calls.append("shortlist")
                start_stage(root, "shortlist")
                finish_stage(root, "shortlist", succeeded=True)
                return {"selected": 10, "total": 100}

            ranking = [
                {"sequestra_candidate_id": f"d{index}", "sequestra_affinity_probability": 1-index/1000, "sequestra_rank": index}
                for index in range(1, 101)
            ]
            with (
                patch("sequestra.orchestrator.run_boltzgen", complete("boltzgen")),
                patch("sequestra.orchestrator.create_shortlist", shortlist),
                patch("sequestra.orchestrator.load_ranked_candidates", return_value=(ranking, "affinity_probability_binary", Path("metrics.csv"))),
                patch("sequestra.orchestrator.run_boltz2_affinity", complete("boltz2_affinity")),
                patch("sequestra.orchestrator.run_structural_confidence", complete("structural_confidence")),
                patch("sequestra.orchestrator.run_reports", complete("reports")),
            ):
                self.assertEqual(
                    run_workflow(root, config_path=config, shortlist_percentage=10, emit=lambda _: None),
                    0,
                )
            self.assertEqual(
                calls,
                ["boltzgen", "shortlist", "boltz2_affinity", "structural_confidence", "reports"],
            )
            self.assertIsNone(next((s for s in load_run_state(root)["stages"] if s["status"] not in {"completed", "skipped"}), None))


if __name__ == "__main__":
    unittest.main()
