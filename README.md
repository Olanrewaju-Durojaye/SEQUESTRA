# SEQUESTRA

<img width="2172" height="724" alt="sequestra-logo" src="https://github.com/user-attachments/assets/83d8eb16-1f9b-4d92-80de-bdf9aed510f9" />

**Reference-benchmarked orchestration for small-molecule-binding protein design**

SEQUESTRA connects BoltzGen, DLKcat, CatPred and Boltz-2 in a resumable,
auditable workflow. It is designed for researchers who want guided operation
without hiding the scientific decision rules or raw model outputs.

SEQUESTRA directly orchestrates BoltzGen, DLKcat, CatPred and Boltz-2 while keeping each scientific tool in its own isolated Conda environment. Users install these scientific tools separately, and SEQUESTRA does not modify their installations.

> Computational predictions are prioritization evidence, not experimental
> proof of binding, catalytic activity, safety or biological function.

## Workflow

1. Generate protein–ligand designs with BoltzGen.
2. Rank designs by BoltzGen affinity probability and retain a user-selected
   top percentage.
3. In catalytic-reference mode, screen catalytic liability with DLKcat and
   CatPred. In non-catalytic-reference mode this screen is optional and
   informational only.
4. Benchmark each survivor against the reference using Boltz-2 binder
   probability.
5. Evaluate structural confidence only for affinity-qualified candidates.
6. Produce consolidated tables, plots, FASTA files, structures, provenance and
   a final report. A zero-survivor result is recorded as a valid outcome.

## Two reference modes

| Mode | Catalytic-liability stage | Affinity benchmark |
|---|---|---|
| `catalytic-reference` | Required; exclude only when predicted kcat is higher **and** predicted Km is lower than the reference | Candidate probability must be strictly greater than reference + margin |
| `non-catalytic-reference` | Disabled by default; optional results never exclude candidates | Candidate probability must be strictly greater than reference + margin |

Structural confidence uses three configurable SEQUESTRA policy defaults:

- `confidence_score >= 0.70`
- `ligand_iptm >= 0.60`
- `complex_iplddt >= 0.70`

All three must pass. These are transparent workflow criteria, not official
Boltz-2 cutoffs.

## Requirements

- Linux workstation
- Python 3.12 and Conda
- NVIDIA GPU recommended (16 GB VRAM minimum configured by default)
- BoltzGen 0.3.x, DLKcat, CatPred and Boltz-2 installed in separate Conda
  environments
- Internet access when Boltz-2's public MSA server is enabled
- At least 40 GiB free disk and 12 GiB available host RAM recommended

Consult each upstream project for its installation and licensing terms.

## Install SEQUESTRA

```bash
conda activate sequestra
cd /path/to/SEQUESTRA
python -m pip install -e .
sequestra setup
sequestra preflight
```

`sequestra setup` creates a computer-specific configuration at
`~/.config/sequestra/sequestra.config.json`. No personal absolute paths are
stored in this repository.

## Beginner quick start

```bash
sequestra launch
```

Answer the guided questions. When project creation finishes, press **Enter** at:

```text
Start the workflow now? [Y/n]:
```

To start later, or resume after any interruption:

```bash
cd /path/to/the/project
sequestra run
```

Check progress at any time:

```bash
sequestra status
```

That is the complete normal user interface. See [TUTORIAL.md](TUTORIAL.md) for
the installation walkthrough, input preparation, both modes, interpretation,
recovery and data deposition. Advanced commands are documented in
[COMMAND_REFERENCE.md](COMMAND_REFERENCE.md).

## Input validation and expression tags

The project creator verifies the chosen PDB chain, bound ligand CCD, complete
FASTA and coordinate-to-sequence mapping. Residues explicitly labelled
`EXPRESSION TAG` in PDB `SEQADV` records are automatically excluded from
canonical FASTA mapping and recorded in
`inputs/reference_input_modifications.json`; the copied source PDB remains
unchanged.

If a structure was manually edited before SEQUESTRA was run, preserve both the
original and edited structures and document every removed residue. In the
non-catalytic 6O0F pilot, the terminal expression-tag sequence `SNSLEVLFQ` was
manually removed from chain A before project creation. That pilot-specific
modification must remain in its Zenodo provenance even though v1.0.0 handles
properly annotated expression tags automatically.

## Principal outputs

| Stage | Main outputs |
|---|---|
| Shortlist | `shortlist/selection/complete_affinity_ranking.csv`, selected structures and FASTA |
| Catalytic liability | `catalytic_liability/results/catalytic_liability.csv` and SVG |
| Boltz-2 affinity | `boltz2_affinity/results/boltz2_affinity.csv`, SVG and qualified FASTA |
| Structural confidence | `structural_confidence/results/structural_confidence.csv`, SVG and qualified FASTA |
| Final reporting | `reports/workflow_summary.md`, final CSV/FASTA and reproducibility manifest |

Raw Boltz-2 structures and confidence JSON files remain under
`boltz2_affinity/raw/<candidate_id>/`. The `source_structure` column in the
structural-confidence CSV gives the exact structure path.

## Reproducibility and safety

- Stage state is written atomically and every expensive stage is resumable.
- Completed model outputs are reused instead of regenerated.
- BoltzGen uses zero DataLoader subprocesses and one analysis process by
  default, with a configurable available-RAM safety floor.
- External tools run sequentially; preflight never launches a scientific model.
- Commands, thresholds, source checksums, outputs and decisions are retained.

## Development

```bash
python -m unittest discover -s tests -v
python -m build
```

## Citation, license and support

Citation metadata are provided in [CITATION.cff](CITATION.cff). SEQUESTRA is
released under the [MIT License](LICENSE). Please report reproducible software
problems through the repository issue tracker without uploading confidential
sequences or unpublished structures.
