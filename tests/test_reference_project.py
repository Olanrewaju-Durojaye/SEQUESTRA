from __future__ import annotations

import csv
import json
from pathlib import Path
import tempfile
import unittest

from sequestra.project import initialize_reference_project
from sequestra.reference import (
    expression_tag_residues, inspect_pdb, map_resolved_to_complete, read_fasta,
)


def atom(serial: int, atom_name: str, residue: str, chain: str, number: int) -> str:
    return (
        f"ATOM  {serial:5d} {atom_name:^4s} {residue:>3s} {chain}{number:4d}    "
        "   0.000   0.000   0.000  1.00 20.00           C"
    )


def hetatm(serial: int, residue: str, chain: str, number: int) -> str:
    return (
        f"HETATM{serial:5d}  C1  {residue:>3s} {chain}{number:4d}    "
        "   0.000   0.000   0.000  1.00 20.00           C"
    )


class ReferenceProjectTests(unittest.TestCase):
    def _inputs(self, root: Path) -> tuple[Path, Path]:
        pdb = root / "reference.pdb"
        pdb.write_text("\n".join([
            atom(1, "CA", "ALA", "A", 1),
            atom(2, "CA", "ASP", "A", 3),
            atom(3, "CA", "CYS", "B", 1),
            hetatm(4, "PMM", "A", 901),
            "END",
        ]) + "\n")
        fasta = root / "complete.fasta"
        fasta.write_text(">test|Chain A\nACD\n")
        return pdb, fasta

    def test_maps_missing_internal_residue(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pdb, _ = self._inputs(Path(directory))
            residues, ligands = inspect_pdb(pdb, "A")
            rows = map_resolved_to_complete("ACD", residues)
            self.assertEqual(ligands, {"PMM"})
            self.assertEqual([row["resolved"] for row in rows], [True, False, True])

    def test_excludes_seqadv_expression_tag_from_canonical_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            pdb = Path(directory) / "tagged.pdb"
            seqadv = (
                "SEQADV 6O0F SER A 826  UNP  P31226    SER   845  EXPRESSION TAG"
            )
            pdb.write_text("\n".join([
                seqadv,
                atom(1, "CA", "THR", "A", 824),
                atom(2, "CA", "CYS", "A", 825),
                atom(3, "CA", "SER", "A", 826),
                hetatm(4, "9SL", "A", 901),
                "END",
            ]) + "\n")
            self.assertEqual(expression_tag_residues(pdb, "A"), {("826", "")})
            residues, ligands = inspect_pdb(pdb, "A")
            self.assertEqual("".join(item.amino_acid for item in residues), "TC")
            self.assertEqual(ligands, {"9SL"})
            rows = map_resolved_to_complete("ATC", residues)
            self.assertEqual([row["resolved"] for row in rows], [False, True, True])

            fasta = Path(directory) / "tagged.fasta"
            fasta.write_text(">canonical\nATC\n")
            project = Path(directory) / "tagged-project"
            result = initialize_reference_project(
                project_dir=project,
                project_name="tagged",
                reference_mode="non-catalytic-reference",
                reference_pdb=pdb,
                reference_chain="A",
                complete_fasta=fasta,
                ligand_name="saxitoxin",
                ligand_ccd="9SL",
                ligand_smiles="C",
                ligand_inchikey=None,
                design_min_length=50,
                design_max_length=80,
                number_of_designs=10,
            )
            self.assertEqual(result["excluded_expression_tag_residue_count"], 1)
            provenance = json.loads(
                (project / "inputs/reference_input_modifications.json").read_text()
            )
            self.assertEqual(
                provenance["coordinate_validation_exclusions"][0]["pdb_author_number"],
                "826",
            )
            self.assertIn(seqadv, (project / "inputs/reference_structure.pdb").read_text())

    def test_initializes_self_contained_project(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb, fasta = self._inputs(root)
            project = root / "project"
            result = initialize_reference_project(
                project_dir=project,
                project_name="pilot",
                reference_mode="catalytic-reference",
                reference_pdb=pdb,
                reference_chain="A",
                complete_fasta=fasta,
                ligand_name="test ligand",
                ligand_ccd="PMM",
                ligand_smiles="COP(=O)(O)O",
                ligand_inchikey="TEST",
                design_min_length=50,
                design_max_length=80,
                number_of_designs=100,
            )
            self.assertEqual(result["unresolved_count"], 1)
            specification = json.loads((project / "configuration/project_specification.yaml").read_text())
            self.assertEqual(
                specification["design"]["initial_ranking_metric_family"],
                "boltzgen_affinity_probability",
            )
            self.assertEqual(specification["reference"]["mode"], "catalytic-reference")
            self.assertEqual(specification["decision_rules"]["kinetic_threshold_source"], "reference_predictions")
            self.assertEqual(specification["decision_rules"]["top_fraction_selection_status"], "pending_after_boltzgen_ranking")
            self.assertTrue((project / "configuration/tool_commands.json").is_file())
            self.assertTrue((project / "configuration/run_state.json").is_file())
            self.assertTrue((project / "inputs/ligand_record.yaml").is_file())
            self.assertTrue((project / "structural_confidence").is_dir())
            with (project / "inputs/reference_residue_mapping.csv").open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 3)

    def test_refuses_nonempty_destination(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            pdb, fasta = self._inputs(root)
            project = root / "project"
            project.mkdir()
            (project / "existing.txt").write_text("preserve me")
            with self.assertRaises(FileExistsError):
                initialize_reference_project(
                    project_dir=project,
                    project_name="pilot",
                    reference_mode="catalytic-reference",
                    reference_pdb=pdb,
                    reference_chain="A",
                    complete_fasta=fasta,
                    ligand_name="test ligand",
                    ligand_ccd="PMM",
                    ligand_smiles="COP(=O)(O)O",
                    ligand_inchikey=None,
                    design_min_length=50,
                    design_max_length=80,
                    number_of_designs=100,
                )
            self.assertEqual((project / "existing.txt").read_text(), "preserve me")


if __name__ == "__main__":
    unittest.main()
