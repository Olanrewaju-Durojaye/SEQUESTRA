"""Reference structure and sequence validation."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True, slots=True)
class FastaRecord:
    header: str
    sequence: str


@dataclass(frozen=True, slots=True)
class ResolvedResidue:
    chain: str
    author_number: str
    insertion_code: str
    residue_name: str
    amino_acid: str


def expression_tag_residues(path: Path, chain: str) -> set[tuple[str, str]]:
    """Return author residue identifiers explicitly annotated as expression tags.

    PDB ``SEQADV`` records distinguish engineered expression-tag residues from
    the canonical biological sequence.  These residues may possess coordinates,
    but they must not make an otherwise correct UniProt FASTA fail validation.
    """
    tagged: set[tuple[str, str]] = set()
    for line in path.expanduser().read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("SEQADV") or "EXPRESSION TAG" not in line.upper():
            continue
        line_chain = line[16:17] if len(line) >= 23 else ""
        author_number = line[18:22].strip() if len(line) >= 23 else ""
        insertion_code = line[22:23].strip() if len(line) >= 23 else ""
        if not line_chain or not author_number:
            fields = line.split()
            if len(fields) >= 5:
                line_chain, author_number, insertion_code = fields[3], fields[4], ""
        if line_chain == chain and author_number:
            tagged.add((author_number, insertion_code))
    return tagged


def read_fasta(path: Path) -> list[FastaRecord]:
    records: list[FastaRecord] = []
    header: str | None = None
    sequence: list[str] = []
    for raw_line in path.expanduser().read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append(_validated_record(header, "".join(sequence)))
            header = line[1:].strip()
            sequence = []
        else:
            if header is None:
                raise ValueError(f"FASTA sequence appears before a header in {path}")
            sequence.append(line.replace(" ", "").upper())
    if header is not None:
        records.append(_validated_record(header, "".join(sequence)))
    if not records:
        raise ValueError(f"no FASTA records found in {path}")
    return records


def _validated_record(header: str, sequence: str) -> FastaRecord:
    if not sequence:
        raise ValueError(f"FASTA record {header!r} has no sequence")
    invalid = sorted(set(sequence) - STANDARD_AA)
    if invalid:
        raise ValueError(f"FASTA record {header!r} contains unsupported residues: {invalid}")
    return FastaRecord(header=header, sequence=sequence)


def select_complete_sequence(records: list[FastaRecord], chain: str) -> FastaRecord:
    if len(records) == 1:
        return records[0]
    chain_token = f"Chain {chain}"
    chains_token = f"Chains {chain}"
    matches = [r for r in records if chain_token in r.header or chains_token in r.header]
    if len(matches) == 1:
        return matches[0]
    raise ValueError(
        f"could not select one complete FASTA record for chain {chain!r}; "
        f"found {len(matches)} matches among {len(records)} records"
    )


def inspect_pdb(path: Path, chain: str) -> tuple[list[ResolvedResidue], set[str]]:
    residues: list[ResolvedResidue] = []
    seen: set[tuple[str, str, str, str]] = set()
    ligands: set[str] = set()
    excluded_expression_tags = expression_tag_residues(path, chain)
    for line in path.expanduser().read_text(encoding="utf-8", errors="replace").splitlines():
        record = line[:6].strip()
        if record not in {"ATOM", "HETATM"} or len(line) < 27:
            continue
        line_chain = line[21:22]
        if line_chain != chain:
            continue
        residue_name = line[17:20].strip().upper()
        author_number = line[22:26].strip()
        insertion_code = line[26:27].strip()
        if record == "ATOM" and (author_number, insertion_code) in excluded_expression_tags:
            continue
        if record == "HETATM":
            if residue_name not in {"HOH", "WAT", "DOD"}:
                ligands.add(residue_name)
            continue
        key = (line_chain, author_number, insertion_code, residue_name)
        if key in seen:
            continue
        seen.add(key)
        amino_acid = AA3_TO_1.get(residue_name)
        if amino_acid is None:
            raise ValueError(
                f"unsupported polymer residue {residue_name} at chain {chain} "
                f"position {author_number}{insertion_code}"
            )
        residues.append(
            ResolvedResidue(
                chain=line_chain,
                author_number=author_number,
                insertion_code=insertion_code,
                residue_name=residue_name,
                amino_acid=amino_acid,
            )
        )
    if not residues:
        raise ValueError(f"no resolved polymer residues found for chain {chain!r} in {path}")
    return residues, ligands


def map_resolved_to_complete(
    complete_sequence: str,
    residues: list[ResolvedResidue],
) -> list[dict[str, str | int | bool]]:
    resolved_sequence = "".join(residue.amino_acid for residue in residues)
    mapping: list[int] = []
    complete_index = 0
    for resolved_index, amino_acid in enumerate(resolved_sequence, 1):
        found = complete_sequence.find(amino_acid, complete_index)
        if found < 0:
            context = resolved_sequence[max(0, resolved_index - 6):resolved_index + 5]
            raise ValueError(
                "the coordinate-derived sequence is not a subsequence of the complete "
                f"sequence near resolved residue {resolved_index} ({context})"
            )
        mapping.append(found)
        complete_index = found + 1

    by_complete = {complete_position: residues[i] for i, complete_position in enumerate(mapping)}
    rows: list[dict[str, str | int | bool]] = []
    for zero_based, amino_acid in enumerate(complete_sequence):
        residue = by_complete.get(zero_based)
        rows.append(
            {
                "complete_position": zero_based + 1,
                "amino_acid": amino_acid,
                "resolved": residue is not None,
                "pdb_chain": residue.chain if residue else "",
                "pdb_author_number": residue.author_number if residue else "",
                "pdb_insertion_code": residue.insertion_code if residue else "",
                "pdb_residue_name": residue.residue_name if residue else "",
            }
        )
    return rows
