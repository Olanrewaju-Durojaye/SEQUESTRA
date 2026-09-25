"""SEQUESTRA command-line entry point."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import shlex

from .preflight import print_report, run_preflight, write_json_report
from .project import initialize_reference_project
from .decision_rules import PROJECT_MODES
from .wizard import collect_answers, create_project, load_answers
from .run_state import (
    initialize_from_project, load_run_state, next_resumable_stage,
    recover_interrupted, status_lines,
)
from .boltzgen_runner import build_boltzgen_command, prepare_design_specification, run_boltzgen
from .shortlist import create_shortlist, load_ranked_candidates, recommendation_for_total, selection_count
from .catalytic_liability import _legacy, build_predictor_commands, prepare_catalytic_liability, run_catalytic_liability
from .boltz2_affinity import build_boltz2_commands, prepare_boltz2_affinity, run_boltz2_affinity
from .structural_confidence import run_structural_confidence
from .orchestrator import find_config, find_project, run_workflow
from .setup import configure_workstation


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="sequestra",
        description="SEQUESTRA guided project application",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("setup", help="configure this computer once")
    preflight = subparsers.add_parser(
        "preflight",
        help="perform a read-only installation and resource audit",
    )
    preflight.add_argument("--config", type=Path)
    preflight.add_argument("--json-output", type=Path)
    initialize = subparsers.add_parser(
        "init-reference",
        help="validate and initialize a reference-benchmark project",
    )
    initialize.add_argument("--project-dir", type=Path, required=True)
    initialize.add_argument("--project-name", required=True)
    initialize.add_argument(
        "--reference-mode",
        choices=PROJECT_MODES,
        required=True,
    )
    initialize.add_argument("--reference-pdb", type=Path, required=True)
    initialize.add_argument("--reference-chain", required=True)
    initialize.add_argument("--complete-fasta", type=Path, required=True)
    initialize.add_argument("--ligand-name", required=True)
    initialize.add_argument("--ligand-ccd", required=True)
    initialize.add_argument("--ligand-smiles", required=True)
    initialize.add_argument("--ligand-inchikey")
    initialize.add_argument("--design-min-length", type=int, required=True)
    initialize.add_argument("--design-max-length", type=int, required=True)
    initialize.add_argument("--number-of-designs", type=int, required=True)
    initialize.add_argument(
        "--top-fraction", type=float,
        help="optional at initialization; normally selected after BoltzGen ranking",
    )
    initialize.add_argument("--affinity-comparison-margin", type=float, default=0.0)
    initialize.add_argument("--enable-catalytic-liability", action="store_true")
    initialize.add_argument("--kcat-unit")
    initialize.add_argument("--km-unit")
    launch = subparsers.add_parser(
        "launch",
        help="open the guided SEQUESTRA project-creation application",
    )
    launch.add_argument(
        "--answers", type=Path,
        help="JSON answers file for reproducible non-interactive creation",
    )
    launch.add_argument(
        "--yes", action="store_true",
        help="create without the final confirmation (requires --answers)",
    )
    launch.add_argument("--start", action="store_true", help="start the workflow after creation")
    launch.add_argument("--shortlist-percentage", type=float)
    launch.add_argument("--devices", type=int, default=1)
    run = subparsers.add_parser("run", help="start or resume the complete guided workflow")
    run.add_argument("--project", type=Path, help="project directory; defaults to the current directory")
    run.add_argument("--config", type=Path, help="optional configuration override")
    run.add_argument("--devices", type=int, default=1)
    run.add_argument("--shortlist-percentage", type=float)
    workflow_init = subparsers.add_parser(
        "workflow-init",
        help="add resumable run state to an existing SEQUESTRA project",
    )
    workflow_init.add_argument("--project", type=Path, required=True)
    status = subparsers.add_parser("status", help="show checkpoint status for a project")
    status.add_argument("--project", type=Path)
    status.add_argument("--json", action="store_true", help="print machine-readable state")
    resume = subparsers.add_parser(
        "resume",
        help="recover interrupted state and report the next stage to execute",
    )
    resume.add_argument("--project", type=Path)
    resume.add_argument("--config", type=Path)
    resume.add_argument("--devices", type=int, default=1)
    resume.add_argument("--shortlist-percentage", type=float)
    boltzgen_prepare = subparsers.add_parser(
        "boltzgen-prepare",
        help="generate the BoltzGen small-molecule design specification",
    )
    boltzgen_prepare.add_argument("--project", type=Path, required=True)
    boltzgen_run = subparsers.add_parser(
        "boltzgen-run",
        help="execute the checkpointed BoltzGen generation stage",
    )
    boltzgen_run.add_argument("--project", type=Path, required=True)
    boltzgen_run.add_argument("--config", type=Path)
    boltzgen_run.add_argument("--devices", type=int, default=1)
    boltzgen_run.add_argument(
        "--smoke-test", action="store_true",
        help="generate at most two designs and leave the full stage pending",
    )
    boltzgen_run.add_argument(
        "--yes", action="store_true", help="execute without interactive confirmation"
    )
    shortlist = subparsers.add_parser(
        "shortlist",
        help="rank by BoltzGen affinity probability and retain a user-selected percentage",
    )
    shortlist.add_argument("--project", type=Path, required=True)
    shortlist.add_argument("--percentage", type=float)
    shortlist.add_argument("--yes", action="store_true", help="select without final confirmation")
    catalytic = subparsers.add_parser(
        "catalytic-liability",
        help="run direct, independently resumable DLKcat and CatPred prediction",
    )
    catalytic.add_argument("--project", type=Path, required=True)
    catalytic.add_argument(
        "--config", type=Path,
        help="SEQUESTRA workstation configuration",
    )
    catalytic.add_argument(
        "--prepare-only", action="store_true",
        help="validate and prepare both predictor jobs without model inference",
    )
    catalytic.add_argument("--yes", action="store_true", help="execute without confirmation")
    affinity = subparsers.add_parser(
        "boltz2-affinity",
        help="predict candidate and reference affinities directly with Boltz-2",
    )
    affinity.add_argument("--project", type=Path, required=True)
    affinity.add_argument("--config", type=Path)
    affinity.add_argument("--prepare-only", action="store_true", help="write and validate Boltz-2 YAML inputs without inference")
    affinity.add_argument("--yes", action="store_true", help="execute without confirmation")
    confidence = subparsers.add_parser(
        "structural-confidence",
        help="assess affinity-qualified candidates using preserved Boltz-2 confidence outputs",
    )
    confidence.add_argument("--project", type=Path, required=True)
    confidence.add_argument("--config", type=Path)
    confidence.add_argument("--yes", action="store_true", help="evaluate without confirmation")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    if args.command == "setup":
        try:
            configure_workstation()
            return 0
        except (ValueError, FileNotFoundError, OSError, json.JSONDecodeError) as error:
            print(f"Setup failed: {error}")
            return 2
    if args.command == "preflight":
        try:
            config = find_config(args.config)
        except FileNotFoundError as error:
            print(f"Unable to locate configuration: {error}")
            return 2
        report = run_preflight(config)
        print_report(report)
        if args.json_output:
            write_json_report(report, args.json_output)
            print(f"JSON report: {args.json_output.expanduser().resolve()}")
        return 0 if report["overall_status"] != "fail" else 2
    if args.command == "init-reference":
        result = initialize_reference_project(
            project_dir=args.project_dir,
            project_name=args.project_name,
            reference_mode=args.reference_mode,
            reference_pdb=args.reference_pdb,
            reference_chain=args.reference_chain,
            complete_fasta=args.complete_fasta,
            ligand_name=args.ligand_name,
            ligand_ccd=args.ligand_ccd,
            ligand_smiles=args.ligand_smiles,
            ligand_inchikey=args.ligand_inchikey,
            design_min_length=args.design_min_length,
            design_max_length=args.design_max_length,
            number_of_designs=args.number_of_designs,
            top_fraction=args.top_fraction,
            affinity_comparison_margin=args.affinity_comparison_margin,
            catalytic_liability_enabled=args.enable_catalytic_liability,
            kcat_unit=args.kcat_unit,
            km_unit=args.km_unit,
        )
        print("SEQUESTRA reference project initialized")
        print(f"Project: {result['project_dir']}")
        print(f"Complete sequence: {result['complete_length']} residues")
        print(f"Resolved sequence: {result['resolved_length']} residues")
        print(f"Unresolved positions: {result['unresolved_count']}")
        print(f"Bound ligand verified: {result['ligand_ccd']}")
        print(f"Reference mode: {result['reference_mode']}")
        return 0
    if args.command == "launch":
        if args.yes and not args.answers:
            parser.error("--yes requires --answers")
        answers = load_answers(args.answers) if args.answers else collect_answers()
        try:
            result = create_project(answers=answers, assume_yes=args.yes)
        except (ValueError, FileNotFoundError, FileExistsError, json.JSONDecodeError) as error:
            print(f"Project creation failed: {error}")
            return 2
        if result is None:
            return 1
        project = Path(result["project_dir"])
        if args.start:
            return run_workflow(
                project,
                devices=args.devices,
                shortlist_percentage=args.shortlist_percentage,
            )
        if args.yes:
            print(f"To start later: cd {project} && sequestra run")
            return 0
        answer = input("Start the workflow now? [Y/n]: ").strip().lower()
        if answer in {"", "y", "yes"}:
            return run_workflow(project, devices=args.devices)
        print(f"To start later: cd {project} && sequestra run")
        return 0
    if args.command == "run":
        try:
            return run_workflow(
                args.project,
                config_path=args.config,
                devices=args.devices,
                shortlist_percentage=args.shortlist_percentage,
            )
        except (ValueError, FileNotFoundError, KeyError, RuntimeError, OSError, json.JSONDecodeError) as error:
            print(f"Unable to run workflow: {error}")
            return 2
    if args.command == "workflow-init":
        try:
            state = initialize_from_project(args.project)
        except (ValueError, FileNotFoundError, FileExistsError, KeyError) as error:
            print(f"Workflow initialization failed: {error}")
            return 2
        print("Resumable workflow state initialized.")
        print(f"Next resumable stage: {next_resumable_stage(state)}")
        return 0
    if args.command == "status":
        try:
            project = find_project(args.project) if args.project else find_project()
            state = load_run_state(project)
        except (ValueError, FileNotFoundError) as error:
            print(f"Unable to read workflow status: {error}")
            return 2
        if args.json:
            print(json.dumps(state, indent=2))
        else:
            print("\n".join(status_lines(state)))
        return 0
    if args.command == "resume":
        try:
            project = find_project(args.project) if args.project else find_project()
            recover_interrupted(project)
            return run_workflow(
                project,
                config_path=args.config,
                devices=args.devices,
                shortlist_percentage=args.shortlist_percentage,
            )
        except (ValueError, FileNotFoundError, KeyError, RuntimeError, OSError, json.JSONDecodeError) as error:
            print(f"Unable to resume workflow: {error}")
            return 2
    if args.command == "boltzgen-prepare":
        try:
            path = prepare_design_specification(args.project)
        except (ValueError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
            print(f"BoltzGen preparation failed: {error}")
            return 2
        print(f"BoltzGen design specification: {path}")
        return 0
    if args.command == "boltzgen-run":
        if args.devices < 1:
            parser.error("--devices must be at least 1")
        try:
            config = find_config(args.config)
            command = build_boltzgen_command(
                args.project,
                config_path=config,
                smoke_test=args.smoke_test,
                devices=args.devices,
            )
        except (ValueError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
            print(f"BoltzGen command preparation failed: {error}")
            return 2
        print("BoltzGen command preview:")
        print(shlex.join(command))
        if args.smoke_test:
            print("Smoke-test mode: at most two designs; the full stage will remain pending.")
        if not args.yes:
            answer = input("Execute BoltzGen now? [y/N]: ").strip().lower()
            if answer not in {"y", "yes"}:
                print("BoltzGen execution cancelled; no stage was started.")
                return 1
        try:
            return run_boltzgen(
                args.project,
                config_path=config,
                smoke_test=args.smoke_test,
                devices=args.devices,
            )
        except (ValueError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
            print(f"BoltzGen execution failed before launch: {error}")
            return 2
    if args.command == "shortlist":
        try:
            state = load_run_state(args.project)
            expected = next_resumable_stage(state)
            if expected != "shortlist":
                raise ValueError(
                    f"shortlist cannot run yet; next resumable stage is {expected or 'none'}"
                )
            ranked, metric, _ = load_ranked_candidates(args.project)
        except (ValueError, FileNotFoundError, KeyError) as error:
            print(f"Unable to prepare shortlist: {error}")
            return 2
        total = len(ranked)
        recommendation = recommendation_for_total(total)
        print(f"Candidates available: {total}")
        print(f"Ranking metric only: {metric} (higher is better)")
        print("Top candidates:")
        for row in ranked[: min(20, total)]:
            print(
                f"  {int(row['sequestra_rank']):4d}  "
                f"{row['sequestra_candidate_id']}  "
                f"{float(row['sequestra_affinity_probability']):.6f}"
            )
        print(
            "Non-binding recommendation: "
            f"{recommendation['minimum_percent']}-{recommendation['maximum_percent']}%"
        )
        percentage = args.percentage
        if percentage is None:
            while True:
                raw = input("Percentage to carry forward: ").strip()
                try:
                    percentage = float(raw)
                    count = selection_count(total, percentage)
                    break
                except ValueError as error:
                    print(error)
        else:
            try:
                count = selection_count(total, percentage)
            except ValueError as error:
                print(f"Invalid shortlist percentage: {error}")
                return 2
        print(f"Selection: top {percentage:g}% = {count} candidate(s), using ceiling rounding.")
        print("Tie policy: exact ceiling count; equal scores are ordered by candidate ID.")
        if not args.yes:
            answer = input("Create this shortlist? [y/N]: ").strip().lower()
            if answer not in {"y", "yes"}:
                print("Shortlist creation cancelled; no files were written.")
                return 1
        try:
            result = create_shortlist(args.project, percentage)
        except (ValueError, FileNotFoundError, KeyError, OSError) as error:
            print(f"Shortlist creation failed: {error}")
            return 2
        print(
            f"Shortlist completed: {result['selected']} of {result['total']} candidates retained."
        )
        print(f"Ranking metric used: {result['metric']}")
        print(f"Results: {args.project.expanduser().resolve() / 'shortlist/selection'}")
        return 0
    if args.command == "catalytic-liability":
        try:
            config = find_config(args.config)
            prepare_catalytic_liability(args.project)
            legacy = _legacy(args.project.expanduser().resolve())
            commands = {} if legacy else build_predictor_commands(args.project, sequestra_config_path=config)
        except (ValueError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
            print(f"Catalytic-liability preparation failed: {error}")
            return 2
        print("Catalytic-liability command preview:")
        if legacy:
            print(f"migration: reuse completed v0.8.0 predictions from {legacy[0]}")
        else:
            for name, (command, output) in commands.items():
                print(f"{name}: {shlex.join(command)}")
                print(f"  expected output: {output}")
        if args.prepare_only:
            print("Prepare-only mode: DLKcat and CatPred inference will not run.")
        if not args.yes:
            answer = input("Execute this catalytic-liability stage now? [y/N]: ").strip().lower()
            if answer not in {"y", "yes"}:
                print("Catalytic-liability execution cancelled; no stage was started.")
                return 1
        try:
            return run_catalytic_liability(
                args.project, sequestra_config_path=config,
                prepare_only=args.prepare_only,
            )
        except (ValueError, FileNotFoundError, KeyError, RuntimeError, OSError) as error:
            print(f"Catalytic-liability execution failed: {error}")
            return 2
    if args.command == "boltz2-affinity":
        try:
            config = find_config(args.config)
            prepared = prepare_boltz2_affinity(args.project)
            commands = build_boltz2_commands(args.project, config_path=config)
        except (ValueError, FileNotFoundError, KeyError, json.JSONDecodeError) as error:
            print(f"Boltz-2 affinity preparation failed: {error}")
            return 2
        print("Boltz-2 affinity command preview:")
        for identifier in prepared["record_ids"]:
            print(f"{identifier}: {shlex.join(commands[identifier])}")
        if args.prepare_only:
            print("Prepare-only mode: Boltz-2 inference will not run.")
        if not args.yes:
            answer = input("Execute this Boltz-2 affinity stage now? [y/N]: ").strip().lower()
            if answer not in {"y", "yes"}:
                print("Boltz-2 affinity execution cancelled; no stage was started.")
                return 1
        try:
            return run_boltz2_affinity(args.project, config_path=config, prepare_only=args.prepare_only)
        except (ValueError, FileNotFoundError, KeyError, RuntimeError, OSError) as error:
            print(f"Boltz-2 affinity execution failed: {error}")
            return 2
    if args.command == "structural-confidence":
        if not args.yes:
            answer = input("Evaluate preserved Boltz-2 structural confidence now? [y/N]: ").strip().lower()
            if answer not in {"y", "yes"}:
                print("Structural-confidence assessment cancelled; no stage was started.")
                return 1
        try:
            return run_structural_confidence(args.project, config_path=find_config(args.config))
        except (ValueError, FileNotFoundError, KeyError, RuntimeError, OSError) as error:
            print(f"Structural-confidence assessment failed: {error}")
            return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
