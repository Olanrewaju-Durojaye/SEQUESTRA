from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from sequestra.boltzgen_runner import (
    build_boltzgen_command, prepare_design_specification, run_boltzgen,
)
from sequestra.run_state import initialize_run_state, load_run_state


class BoltzGenRunnerTests(unittest.TestCase):
    def _project(self, root: Path, requested: int = 10) -> tuple[Path, Path]:
        project = root / "project"
        (project / "configuration").mkdir(parents=True)
        (project / "boltzgen").mkdir()
        (project / "logs").mkdir()
        specification = {
            "design": {
                "minimum_length": 50,
                "maximum_length": 80,
                "number_of_designs": requested,
            },
            "ligand": {"ccd": "PMM"},
        }
        (project / "configuration/project_specification.yaml").write_text(
            json.dumps(specification)
        )
        initialize_run_state(
            project,
            reference_mode="catalytic-reference",
            catalytic_liability_enabled=True,
        )
        environment = root / "boltzgen-env"
        executable = environment / "bin/boltzgen"
        executable.parent.mkdir(parents=True)
        executable.write_text(
            "#!/bin/sh\n"
            "out=''\n"
            "while [ $# -gt 0 ]; do\n"
            "  if [ \"$1\" = '--output' ]; then shift; out=$1; fi\n"
            "  shift\n"
            "done\n"
            "mkdir -p \"$out/final_ranked_designs\"\n"
            "printf 'design_id,affinity_probability_binary\\nA,0.9\\n' > \"$out/final_ranked_designs/all_designs_metrics.csv\"\n"
        )
        executable.chmod(0o755)
        config = root / "sequestra.config.json"
        config.write_text(json.dumps({"environments": {"boltzgen": str(environment)}}))
        return project, config

    def test_design_specification_uses_project_values(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project, _ = self._project(Path(directory))
            text = prepare_design_specification(project).read_text()
            self.assertIn("sequence: 50..80", text)
            self.assertIn("ccd: PMM", text)

    def test_command_uses_small_molecule_protocol_and_reuse(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project, config = self._project(Path(directory), requested=10)
            command = build_boltzgen_command(
                project, config_path=config, smoke_test=True, devices=1
            )
            self.assertIn("protein-small_molecule", command)
            self.assertIn("--reuse", command)
            self.assertEqual(command[command.index("--num_workers") + 1], "0")
            analysis = command.index("--config")
            self.assertEqual(command[analysis + 1:analysis + 3], ["analysis", "num_processes=1"])
            self.assertEqual(command[command.index("--num_designs") + 1], "2")

    def test_config_can_override_resource_safeguards(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project, config = self._project(Path(directory), requested=10)
            payload = json.loads(config.read_text())
            payload["boltzgen"] = {
                "num_workers": 1,
                "minimum_available_ram_gib": 3,
                "memory_poll_seconds": 0.5,
            }
            config.write_text(json.dumps(payload))
            command = build_boltzgen_command(
                project, config_path=config, smoke_test=False, devices=1
            )
            self.assertEqual(command[command.index("--num_workers") + 1], "1")

    def test_smoke_test_leaves_full_stage_pending(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project, config = self._project(Path(directory), requested=10)
            result = run_boltzgen(
                project, config_path=config, smoke_test=True, emit=lambda _: None
            )
            self.assertEqual(result, 0)
            stage = load_run_state(project)["stages"][0]
            self.assertEqual(stage["status"], "pending")
            self.assertEqual(stage["attempts"][0]["outcome"], "smoke_test_completed")

    def test_full_run_completes_stage(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project, config = self._project(Path(directory), requested=2)
            result = run_boltzgen(
                project, config_path=config, smoke_test=False, emit=lambda _: None
            )
            self.assertEqual(result, 0)
            stage = load_run_state(project)["stages"][0]
            self.assertEqual(stage["status"], "completed")


if __name__ == "__main__":
    unittest.main()
