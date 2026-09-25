# SEQUESTRA command reference

The beginner interface is:

```bash
sequestra launch
cd /path/to/project
sequestra run
sequestra status
```

Use the following only for scripted operation, inspection or recovery.

## Configuration

```bash
sequestra setup
sequestra preflight
sequestra preflight --json-output preflight-report.json
```

Configuration precedence is an explicit `--config`, `SEQUESTRA_CONFIG`, then
`~/.config/sequestra/sequestra.config.json`.

## Non-interactive launch

```bash
cp examples/project_answers.example.json my_answers.json
# Edit my_answers.json, then:
sequestra launch --answers my_answers.json --yes
```

## Explicit run

```bash
sequestra run --project /path/to/project --shortlist-percentage 10
```

## Stage preparation and execution

```bash
sequestra boltzgen-prepare --project /path/to/project
sequestra boltzgen-run --project /path/to/project --smoke-test
sequestra shortlist --project /path/to/project --percentage 10 --yes
sequestra catalytic-liability --project /path/to/project --prepare-only --yes
sequestra boltz2-affinity --project /path/to/project --prepare-only --yes
sequestra structural-confidence --project /path/to/project --yes
```

Remove `--prepare-only` to execute the corresponding model stage. Ordinary
users do not need to pass `--config`; the workstation configuration is
discovered automatically.

## Existing early-version projects

```bash
sequestra workflow-init --project /path/to/project
```

This adds a run-state manifest without recreating completed design inputs.

## Decision rules

- Shortlist: descending BoltzGen affinity probability; ceiling-rounded count.
- Catalytic exclusion: candidate kcat strictly higher **and** Km strictly lower
  than reference.
- Affinity: candidate probability strictly greater than reference + margin.
- Structural confidence: every configured minimum must be satisfied.

Raw model outputs remain authoritative and are never replaced by plot-derived
values.
