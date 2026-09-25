from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from sequestra.preflight import _load_config, run_preflight


class PreflightTests(unittest.TestCase):
    def test_rejects_unknown_schema(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            path.write_text('{"schema_version": 99, "environments": {}}')
            with self.assertRaises(ValueError):
                _load_config(path)

    def test_report_declares_no_model_execution(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            env = root / "env"
            (env / "bin").mkdir(parents=True)
            python = env / "bin/python"
            python.write_text("#!/bin/sh\nexit 1\n")
            python.chmod(0o755)
            boltz = env / "bin/boltz"
            boltz.write_text("#!/bin/sh\nexit 0\n")
            boltz.chmod(0o755)
            conda = root / "conda"
            conda.write_text("#!/bin/sh\nexit 0\n")
            conda.chmod(0o755)
            repository = root / "repo"
            repository.mkdir()
            config = {
                "schema_version": 1,
                "conda_executable": str(conda),
                "environments": {
                    "sequestra": str(env),
                    "boltzgen": str(env),
                    "dlkcat": str(env),
                    "catpred": str(env),
                    "boltz2": str(env),
                },
                "repositories": {
                    "sequestra": str(repository),
                    "boltra_reference": str(repository),
                },
                "predictors": {
                    name: {"repository": str(repository), "command": ["{python}", "predict.py"], "output_csv": "{output_dir}/predictions.csv"}
                    for name in ("dlkcat", "catpred")
                },
                "minimum_resources": {},
            }
            config_path = root / "config.json"
            config_path.write_text(json.dumps(config))
            with patch.dict(os.environ, {"CONDA_PREFIX": str(env)}), patch(
                "sequestra.preflight._gpu_info",
                return_value=("warn", None, "test GPU unavailable"),
            ):
                report = run_preflight(config_path)
            self.assertTrue(report["read_only"])
            self.assertFalse(report["scientific_models_executed"])
            self.assertIn(report["overall_status"], {"pass", "warn"})


if __name__ == "__main__":
    unittest.main()
