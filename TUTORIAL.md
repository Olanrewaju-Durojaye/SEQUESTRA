# SEQUESTRA beginner-friendly tutorial

This tutorial assumes no command-line programming experience. Commands shown in
code boxes can be copied into a Linux terminal. Replace text beginning with
`/path/to/` with the location on your own computer.

## 1. What SEQUESTRA does

SEQUESTRA manages a sequence of independent scientific programs. It creates
inputs, launches one stage at a time, checks expected outputs, records decisions
and resumes safely after interruption. It does not claim that a computational
hit is experimentally validated.

Two workflows are available:

- **Catalytic reference:** the natural reference catalyses a reaction involving
  the ligand. Candidate designs with a stronger adverse two-part kinetic signal
  are excluded before affinity benchmarking.
- **Non-catalytic binding reference:** the reference binds the ligand without
  being used as a catalytic benchmark. Catalytic-liability prediction is
  normally disabled.

## 2. Prepare the computer once

Install Miniconda or Anaconda, the NVIDIA driver/CUDA stack required by the
upstream tools, and these isolated environments:

| Environment name | Required software |
|---|---|
| `sequestra` | SEQUESTRA and Python 3.12 |
| `boltzgen` | BoltzGen |
| `dlkcat` | DLKcat and its model resources |
| `catpred` | CatPred, checkpoints and data directory |
| `boltz2` | Boltz-2 (`boltz` executable) |

SEQUESTRA deliberately does not merge these environments. This prevents one
tool's dependencies from silently breaking another tool.

Install this repository:

```bash
conda activate sequestra
cd /path/to/SEQUESTRA
python -m pip install -e .
python -c "import sequestra; print(sequestra.__version__)"
```

The final line should print `1.0.0`.

Configure local paths and audit without running a scientific model:

```bash
sequestra setup
sequestra preflight
```

SEQUESTRA detects Conda environments by name and asks only for directories it
cannot find. Do not begin a large campaign unless the overall status is
`PASS`. A public MSA server may still require reliable internet during Boltz-2.

## 3. Prepare the biological inputs

Create one preparation folder containing:

1. A PDB file containing the selected protein chain and bound ligand.
2. A complete canonical FASTA sequence for that chain.
3. The ligand name, PDB Chemical Component Dictionary (CCD) code, isomeric
   SMILES and, when available, InChIKey.

### PDB checklist

- Confirm the intended chain identifier, for example `A`.
- Confirm that the ligand CCD occurs as a bound non-water component.
- Keep crystallographic cofactors or ions only when scientifically intended.
- Do not silently edit missing residues into an experimental structure.
- Preserve the original downloaded PDB even when preprocessing is required.

### Expression tags and engineered residues

Version 1.0.0 recognizes residues explicitly annotated as `EXPRESSION TAG` in
PDB `SEQADV` records. They are excluded from canonical FASTA mapping and listed
in `inputs/reference_input_modifications.json`. The project retains the source
PDB unchanged.

For the deposited non-catalytic 6O0F pilot, preprocessing was performed before
this automatic behavior was available: terminal chain-A residues
`SNSLEVLFQ` were manually removed. Deposit both the original and edited PDB,
identify the removed residues, give the reason (engineered expression tag), and
state that no canonical biological residues were intentionally removed.

If an unannotated engineered segment causes a sequence mismatch, do not append
it to the canonical FASTA merely to force validation. Verify the construct from
the PDB record or publication, preserve the original, make an auditable copy,
and document the exact edit.

### FASTA checklist

- Use the complete canonical sequence, not only residues visible in the PDB.
- Use the correct isoform and chain.
- Standard one-letter amino-acid characters are required.
- With a multi-record FASTA, ensure the headers identify the relevant chain.

## 4. Create a project

Activate SEQUESTRA and launch the guide:

```bash
conda activate sequestra
sequestra launch
```

Choose `1` for catalytic reference or `2` for non-catalytic binding reference.
The guide asks for the project destination, input files, ligand description,
design lengths, design count and Boltz-2 comparison margin.

For a non-catalytic reference, answer `n` to the optional catalytic-liability
screen unless prediction-only kinetic information is specifically wanted.
Enabling it never changes which candidates advance in that mode.

Review the displayed summary carefully and confirm project creation. The
destination must be new or empty; SEQUESTRA refuses to overwrite an existing
project.

When asked:

```text
Start the workflow now? [Y/n]:
```

press **Enter** to start. To start later, answer `n`.

## 5. Start or resume with one command

Move into the project and run:

```bash
cd /path/to/project
sequestra run
```

You do not need to supply the project path, configuration path, stage name or
GPU number for normal use. Running the same command after a shutdown,
interruption or recoverable network failure resumes from the first incomplete
stage.

Never launch two `sequestra run` commands for the same project simultaneously.

## 6. Choose the shortlist percentage

After BoltzGen finishes, SEQUESTRA reports the number of ranked designs and
asks what percentage to carry forward. With 100 designs, 10% retains 10 and
20% retains 20 candidates.

Selection uses descending BoltzGen `affinity_probability_binary` and ceiling
rounding. The recommendation is guidance, not a scientific requirement. The
confirmed percentage is recorded with the ranking and source checksum.

## 7. Understand each decision stage

### 7.1 Catalytic liability

For a catalytic reference, a candidate is excluded only when both are true:

```text
candidate predicted kcat > reference predicted kcat
AND
candidate predicted Km < reference predicted Km
```

Equality does not trigger an adverse condition. One adverse condition is a
partial-liability signal and the candidate continues. DLKcat and CatPred are
predictions; their output does not establish catalytic function.

### 7.2 Boltz-2 affinity

The reference and each eligible candidate are modelled with the same ligand.
A candidate advances only when:

```text
candidate affinity_probability_binary
>
reference affinity_probability_binary + predeclared margin
```

Equality fails the strict comparison. `affinity_pred_value` is retained but is
not the gate used by SEQUESTRA.

### 7.3 Structural confidence

Boltz-2 is not run again. SEQUESTRA reuses the preserved `model_0` structure
and confidence JSON. The defaults require:

```text
confidence_score >= 0.70
ligand_iptm >= 0.60
complex_iplddt >= 0.70
```

All three must pass. These configurable thresholds are SEQUESTRA screening
policy, not experimentally calibrated or official Boltz-2 boundaries.

If no candidate survives affinity, structural confidence is marked `SKIPPED`.
If candidates enter structural confidence but none pass, the stage is
`COMPLETED` with zero final candidates. Both are valid negative outcomes.

## 8. Monitor the run

From anywhere inside the project:

```bash
sequestra status
```

| State | Meaning |
|---|---|
| `PENDING` | Not started |
| `RUNNING` | Currently running, or interrupted before recovery |
| `COMPLETED` | Finished and validated |
| `FAILED` | Stopped because an error requires attention |
| `SKIPPED` | Intentionally unnecessary for this workflow/outcome |

Logs are stored under `logs/`. Do not delete raw stage directories while a run
is resumable.

## 9. Common warnings and recovery

### ResourceWarning about unclosed files

Older BoltzGen versions may print Python `ResourceWarning` messages while
writing intermediate CIF files. These warnings do not alone mean the model has
failed. Judge completion from the stage exit status and expected outputs.

### Available RAM below the safety floor

SEQUESTRA stops BoltzGen before the whole workstation reaches an out-of-memory
condition. Close other applications or increase swap, then run
`sequestra run`. Completed BoltzGen artifacts are preserved through `--reuse`.

### MSA server timeout or lost internet

Boltz-2 may use `https://api.colabfold.com`. Restore internet access and rerun
`sequestra run`. Completed per-candidate checkpoints are reused.

### DLKcat `GLIBCXX` import error

Test the environment directly:

```bash
conda activate dlkcat
python -c "from rdkit import Chem; print('RDKit import: PASS')"
```

Resolve the Conda runtime-library inconsistency inside the DLKcat environment;
do not replace system libraries or modify unrelated environments.

### Project not found

Enter the project directory (or one of its subdirectories) and run
`sequestra run` again.

## 10. Find and interpret outputs

```text
project/
├── inputs/                       validated input copies and provenance
├── configuration/                frozen specification and run state
├── boltzgen/                      generation outputs
├── shortlist/selection/           ranking and selected designs
├── catalytic_liability/results/   kinetic decisions and plot
├── boltz2_affinity/results/        affinity decisions and plot
├── boltz2_affinity/raw/            Boltz-2 structures/confidence files
├── structural_confidence/results/ confidence decisions and plot
├── reports/                        final or terminal summary
└── logs/                           execution logs
```

Plots use color for scalable classification; candidate names are omitted from
visible labels to avoid crowding. Hover metadata retains identifiers, while the
CSV is the authoritative table.

Locate a candidate's preserved Boltz-2 model with:

```bash
find boltz2_affinity/raw/CANDIDATE_ID \
  -type f \( -name '*model_0.cif' -o -name '*model_0.pdb' \)
```

The `source_structure` column of `structural_confidence.csv` gives the same
exact relative path.

## 11. Prepare a reproducible Zenodo deposition

Deposit enough information to reproduce interpretation without claiming that
large upstream model checkpoints are your own work:

- SEQUESTRA source release and version.
- Project specification, run state, logs and reproducibility manifest.
- Original reference structure and canonical FASTA.
- Any edited reference structure plus a description of every modification.
- Ligand identity fields and input-source citations.
- Consolidated CSV/SVG/FASTA/JSON reports.
- Candidate structures permitted by upstream licenses.
- Hardware, OS, GPU driver and Conda environment exports.
- Exact versions/commits of BoltzGen, DLKcat, CatPred and Boltz-2.
- Dates of public MSA-server access and disclosure of transient failures.
- A statement that predictions lack experimental validation.

For the 6O0F pilot, explicitly report manual removal of chain-A terminal
expression-tag residues `SNSLEVLFQ`, retain both PDB versions, and provide their
SHA-256 checksums. Do not describe this as removal of native saxiphilin
sequence.

## 12. Archive the completed project

After `sequestra status` reports no resumable stage, create a preservation copy
without deleting the working project:

```bash
cd /path/to/project/..
tar -czf project-name.tar.gz project-name/
sha256sum project-name.tar.gz > project-name.tar.gz.sha256
```

Keep the archive checksum with the Zenodo description.

## 13. Advanced and scripted operation

Normal users should prefer `sequestra launch` and `sequestra run`. For
automation and stage-specific commands, see [COMMAND_REFERENCE.md](COMMAND_REFERENCE.md)
and `examples/project_answers.example.json`.

---

## For questions
Please contact the corresponding authors:

Corresponding authors:
    Olanrewaju Ayodeji Durojaye;
    Rachid Daoud

Institution:
    Chemical and Biochemical Sciences, Green Process Engineering,
    University Mohammed VI Polytechnic,
    43150 Ben Guerir, Morocco

Email:
    olanrewaju.ayodeji-durojaye-ext@um6p.ma;
    rachid.daoud@um6p.ma
