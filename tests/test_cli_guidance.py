from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from sequestra.cli import main


class CliGuidanceTests(unittest.TestCase):
    def test_interactive_launch_offers_and_starts_workflow(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            project = Path(directory) / "project"
            with (
                patch("sequestra.cli.collect_answers", return_value={}),
                patch("sequestra.cli.create_project", return_value={"project_dir": str(project)}),
                patch("builtins.input", return_value="y"),
                patch("sequestra.cli.run_workflow", return_value=0) as runner,
            ):
                self.assertEqual(main(["launch"]), 0)
            runner.assert_called_once_with(project, devices=1)


if __name__ == "__main__":
    unittest.main()
