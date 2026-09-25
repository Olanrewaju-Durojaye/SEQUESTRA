from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from sequestra.setup import configure_workstation


class SetupTests(unittest.TestCase):
    def test_setup_discovers_named_environments_and_common_tool_folders(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            home = Path(directory)
            conda = home / "miniconda3/bin/conda"
            conda.parent.mkdir(parents=True)
            conda.write_text("executable")
            envs = []
            for name in ("sequestra", "boltzgen", "dlkcat", "catpred", "boltz2"):
                path = home / f"miniconda3/envs/{name}"
                path.mkdir(parents=True)
                envs.append(str(path))
            (home / "Applications/DLKcat").mkdir(parents=True)
            (home / "Applications/catpred_pipeline/CatPred").mkdir(parents=True)
            (home / "Applications/catpred_pipeline/data").mkdir(parents=True)
            completed = SimpleNamespace(returncode=0, stdout=json.dumps({"envs": envs}))
            with (
                patch("sequestra.setup.Path.home", return_value=home),
                patch("sequestra.setup.shutil.which", return_value=str(conda)),
                patch("sequestra.setup.subprocess.run", return_value=completed),
            ):
                destination = configure_workstation(ask=lambda _: self.fail("unexpected prompt"), emit=lambda _: None)
            config = json.loads(destination.read_text())
            self.assertEqual(config["environments"]["boltz2"], str(home / "miniconda3/envs/boltz2"))
            self.assertEqual(config["predictors"]["dlkcat"]["repository"], str(home / "Applications/DLKcat"))
            self.assertEqual(destination, home / ".config/sequestra/sequestra.config.json")


if __name__ == "__main__":
    unittest.main()
