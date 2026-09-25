from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from sequestra.run_state import (
    finish_stage,
    initialize_run_state,
    load_run_state,
    next_resumable_stage,
    recover_interrupted,
    skip_stage,
    start_stage,
)


class RunStateTests(unittest.TestCase):
    def test_disabled_liability_stage_is_skipped(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            state = initialize_run_state(
                root,
                reference_mode="non-catalytic-reference",
                catalytic_liability_enabled=False,
            )
            liability = next(row for row in state["stages"] if row["name"] == "catalytic_liability")
            self.assertEqual(liability["status"], "skipped")
            self.assertEqual(next_resumable_stage(state), "boltzgen")

    def test_completed_stage_is_skipped_on_resume(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            initialize_run_state(
                root, reference_mode="catalytic-reference", catalytic_liability_enabled=True
            )
            start_stage(root, "boltzgen", command=["boltzgen", "run"])
            state = finish_stage(root, "boltzgen", succeeded=True, outputs=["boltzgen/results.csv"])
            self.assertEqual(next_resumable_stage(state), "shortlist")
            with self.assertRaises(ValueError):
                start_stage(root, "boltzgen")

    def test_running_stage_becomes_interrupted_and_retries(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            initialize_run_state(
                root, reference_mode="catalytic-reference", catalytic_liability_enabled=True
            )
            start_stage(root, "boltzgen")
            recovered = recover_interrupted(root)
            first = recovered["stages"][0]
            self.assertEqual(first["status"], "interrupted")
            self.assertEqual(first["attempts"][0]["outcome"], "interrupted")
            retried = start_stage(root, "boltzgen")
            self.assertEqual(len(retried["stages"][0]["attempts"]), 2)

    def test_failed_stage_is_resumable(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            initialize_run_state(
                root, reference_mode="catalytic-reference", catalytic_liability_enabled=True
            )
            start_stage(root, "boltzgen")
            failed = finish_stage(root, "boltzgen", succeeded=False, message="test failure")
            self.assertEqual(next_resumable_stage(failed), "boltzgen")

    def test_corrupt_state_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            path = root / "configuration/run_state.json"
            path.parent.mkdir()
            path.write_text("not-json")
            with self.assertRaises(ValueError):
                load_run_state(root)

    def test_next_stage_can_be_skipped_with_reason(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            initialize_run_state(root, reference_mode="catalytic-reference", catalytic_liability_enabled=True)
            for name in ("boltzgen", "shortlist", "catalytic_liability", "boltz2_affinity"):
                start_stage(root, name)
                finish_stage(root, name, succeeded=True)
            state = skip_stage(root, "structural_confidence", reason="no_affinity_qualified_candidates")
            structural = next(row for row in state["stages"] if row["name"] == "structural_confidence")
            self.assertEqual(structural["status"], "skipped")
            self.assertEqual(next_resumable_stage(state), "reports")


if __name__ == "__main__":
    unittest.main()
