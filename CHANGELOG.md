# Changelog

## 1.0.0

- Add the beginner-facing `README.md`, comprehensive `TUTORIAL.md`, command
  reference, MIT license and citation metadata for public release.
- Automatically exclude PDB residues explicitly annotated as `EXPRESSION TAG`
  in `SEQADV` records from canonical FASTA mapping.
- Preserve the source PDB unchanged and record mapping exclusions in
  `inputs/reference_input_modifications.json`.
- Remove the developer workstation configuration from the distributed source;
  each computer now creates its own configuration with `sequestra setup`.
- Document the exact affinity and structural-confidence policy, negative-result
  handling, raw 3D model locations, recovery procedures and Zenodo provenance.
- Record the pilot-specific manual removal of `SNSLEVLFQ` from 6O0F chain A as
  a required data-deposition disclosure.

## 0.11.3

- Support CatPred's native `log10km_mean` output schema.
- Convert native CatPred Km predictions from log10(mM) into the Km unit frozen
  in the project (`M`, `mM`, or `uM`).
- Convert MVE variance to a log10-space standard deviation and preserve its
  scale explicitly in the normalized prediction table.
- Cross-check identities using exact sequence and unique `pdbpath` values.
- Reuse existing raw CatPred predictions after post-processing failure instead
  of repeating ten-model inference.

## 0.11.2

- Limit BoltzGen structural analysis to one CPU process instead of its default
  32-process multiprocessing pool, preventing each worker from duplicating the
  100-design analysis dataset in host RAM.
- Retain the v0.11.1 DataLoader and available-RAM safeguards.
- Resume existing campaigns with native `--reuse`; completed generation,
  inverse-folding, folding, design-folding, and affinity outputs remain intact.

## 0.11.1

- Prevent BoltzGen post-affinity analysis from multiplying host-memory use by
  defaulting its PyTorch DataLoader to zero subprocess workers.
- Add a host-RAM watchdog that stops the complete BoltzGen process group before
  Linux reaches global out-of-memory conditions.
- Preserve all partial BoltzGen artifacts for native `--reuse` recovery.
- Recover stale `RUNNING` checkpoints automatically through `sequestra run`.

## 0.11.0

- Added `sequestra run` as the single beginner-facing command for starting and resuming complete workflows.
- Project discovery now works from the project root or any subdirectory.
- Configuration is discovered automatically; ordinary users no longer supply `--config` or `--devices`.
- Added `sequestra setup` for one-time, computer-specific Conda and scientific-tool discovery.
- Interactive project creation now offers to start the workflow immediately.
- The runner pauses after BoltzGen only to collect and confirm the user's shortlist percentage.
- Automatically selects the catalytic or non-catalytic workflow branch and skips completed checkpoints.
- Added complete survivor-path reporting with final CSV/FASTA, Markdown/JSON summary, and checksummed reproducibility manifest.
- `sequestra status` and `sequestra resume` can operate from inside a project without a path.

## 0.10.0

- Added structural-confidence assessment for affinity-qualified candidates only.
- Reuses preserved, highest-ranked Boltz-2 `model_0` structures and confidence JSON files without new inference.
- Added configurable minimum thresholds for overall confidence, protein-ligand ipTM, and interface-weighted pLDDT.
- Produces a consolidated CSV, qualified FASTA, three-panel confidence SVG, and JSON summary.
- Added independently reusable per-candidate confidence checkpoints.
- Supports catalytic- and non-catalytic-reference projects identically after affinity qualification.
- Automatically skips structural confidence when affinity has zero survivors and generates terminal Markdown and JSON reports.
- Records zero survivors as a valid computational negative result rather than a workflow failure.

## 0.9.1

- Added a reference-relative Boltz-2 affinity SVG with numeric 0–1 axis values.
- Shows the strict reference-plus-margin decision threshold and advancement region.
- Uses color-only point classification to remain readable with large candidate sets.
- Keeps candidate IDs, probability values, and decisions in SVG hover metadata and the consolidated CSV.

## 0.9.0

- Added direct Boltz-2 protein–ligand complex and affinity prediction.
- Supports both catalytic- and non-catalytic-reference modes.
- Generates one auditable Boltz-2 YAML input for the reference and each eligible candidate.
- Uses `affinity_probability_binary` with a strict reference-relative comparison; equality does not pass.
- Retains `affinity_pred_value` as a secondary reported value rather than using it for the hit-discovery gate.
- Added independent per-record checkpoints and safe resume behavior.
- Preserves raw Boltz-2 structures, confidence outputs, and affinity JSON files.
- Produces a consolidated decision table and affinity-qualified FASTA for structural-confidence assessment.
- Added direct Boltz executable validation to preflight.

## 0.8.4

- Added SEQUESTRA-owned native DLKcat and CatPred execution adapters.
- Generates each tool's required native input and protein-record files directly.
- Runs predictors in their existing isolated Conda environments without ENZTRA.
- Preserves native raw outputs and writes normalized prediction tables keyed by candidate ID.
- Reconciles outputs by exact protein sequence and rejects missing or ambiguous identities.
- Validates real scripts, data roots, and model resources during preflight.
- Retains independent DLKcat, CatPred, and combination checkpoints and legacy recovery.

## 0.8.3

- Added numeric tick values and light grid lines to both log-scale plot axes.
- Removed candidate-ID text labels from plotted points to prevent crowding.
- Retained candidate identity, values, and decision as SVG hover metadata.

## 0.8.2

- Recognizes the actual ENZTRA `kinetics.csv` columns `kcat_s`, `km_mm`, and `km_sd_total`.
- Prioritizes a combined table containing an explicit `role=reference` row over candidate-only selection tables.
- Completes existing v0.8.0 projects entirely from preserved predictions without launching any predictor.

## 0.8.1

- Removed ENZTRA as a runtime and preflight dependency.
- Added direct configurable DLKcat and CatPred command adapters.
- Added independent DLKcat, CatPred, and combination checkpoints.
- Added identity validation before combining predictor outputs.
- Added transparent recovery from preserved v0.8.0 `enztra_job` predictions.
- Added corrected reversed-rule classification and an exclusion-oriented SVG.
- Preserved compatibility with existing projects and upstream checkpoints.
