from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

from sequestra.decision_rules import shortlist_recommendation
from sequestra.wizard import _new_project_directory, collect_answers, create_project


def atom(serial: int, residue: str, chain: str, number: int) -> str:
    return (
        f"ATOM  {serial:5d}  CA  {residue:>3s} {chain}{number:4d}    "
        "   0.000   0.000   0.000  1.00 20.00           C"
    )


class WizardTests(unittest.TestCase):
    def _answers(self, root: Path) -> dict[str, object]:
        pdb = root / "reference.pdb"
        pdb.write_text("\n".join([
            atom(1, "ALA", "A", 1),
            atom(2, "CYS", "A", 2),
            "HETATM    3  C1  LIG A 901       0.000   0.000   0.000  1.00 20.00           C",
            "END",
        ]) + "\n")
        fasta = root / "reference.fasta"
        fasta.write_text(">reference|Chain A\nAC\n")
        return {
            "project_dir": root / "project",
            "project_name": "guided-test",
            "reference_mode": "non-catalytic-reference",
            "reference_pdb": pdb,
            "reference_chain": "A",
            "complete_fasta": fasta,
            "ligand_name": "test ligand",
            "ligand_ccd": "LIG",
            "ligand_smiles": "CCO",
            "ligand_inchikey": None,
            "design_min_length": 40,
            "design_max_length": 60,
            "number_of_designs": 2000,
            "catalytic_liability_enabled": False,
        }

    def test_scale_aware_recommendations(self) -> None:
        self.assertEqual(shortlist_recommendation(100)["minimum_percent"], 10)
        self.assertEqual(shortlist_recommendation(1000)["minimum_percent"], 5)
        self.assertEqual(shortlist_recommendation(1001)["minimum_percent"], 1)

    def test_guided_creation_keeps_shortlist_pending(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output: list[str] = []
            result = create_project(
                answers=self._answers(root), assume_yes=True, emit=output.append
            )
            self.assertIsNotNone(result)
            self.assertTrue((root / "project/configuration/project_specification.yaml").is_file())
            self.assertEqual(result["shortlist_recommendation"]["minimum_percent"], 1)
            self.assertTrue(any("remains pending" in line for line in output))
            self.assertTrue(any("Non-binding" in line for line in output))

    def test_cancellation_writes_nothing(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            answers = self._answers(root)
            result = create_project(answers=answers, ask=lambda _: "no")
            self.assertIsNone(result)
            self.assertFalse((root / "project").exists())

    def test_nonempty_destination_is_rejected_early(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "input.pdb").write_text("input")
            choices = iter((str(root), str(root / "new-project")))
            selected = _new_project_directory(lambda _: next(choices))
            self.assertEqual(selected, str(root / "new-project"))

    def test_noncatalytic_mode_without_screen_skips_unit_questions(self) -> None:
        responses = iter((
            "2", "/tmp/new-sequestra-project", "pilot", "/tmp/reference.pdb", "A",
            "/tmp/reference.fasta", "ligand", "LIG", "CCO", "", "40", "60",
            "100", "", "n",
        ))
        answers = collect_answers(lambda _: next(responses))
        self.assertFalse(answers["catalytic_liability_enabled"])
        self.assertIsNone(answers["kcat_unit"])
        self.assertIsNone(answers["km_unit"])
        with self.assertRaises(StopIteration):
            next(responses)


if __name__ == "__main__":
    unittest.main()
