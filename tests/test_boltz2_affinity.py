from __future__ import annotations
import csv,json,tempfile,unittest,xml.etree.ElementTree as ET
from pathlib import Path
from sequestra.boltz2_affinity import prepare_boltz2_affinity,run_boltz2_affinity
from sequestra.run_state import initialize_run_state,load_run_state,save_run_state

class Boltz2AffinityTests(unittest.TestCase):
    def project(self,base:Path,mode:str="catalytic-reference")->Path:
        root=base/"project";(root/"configuration").mkdir(parents=True);(root/"inputs").mkdir();(root/"catalytic_liability/results").mkdir(parents=True)
        (root/"inputs/reference_complete.fasta").write_text(">ref\nACDEFG\n")
        (root/"configuration/project_specification.yaml").write_text(json.dumps({"reference":{"mode":mode},"ligand":{"name":"ligand","smiles":"CCO","ccd":"LIG"},"decision_rules":{"affinity_comparison_margin":0.0,"catalytic_liability_enabled":mode=="catalytic-reference"}}))
        (root/"catalytic_liability/results/automatic_survivors.fasta").write_text(">design_a\nHIKLMN\n>design_b\nPQRSTV\n")
        initialize_run_state(root,reference_mode=mode,catalytic_liability_enabled=mode=="catalytic-reference");state=load_run_state(root)
        for stage in state["stages"]:
            if stage["name"] in {"boltzgen","shortlist","catalytic_liability"}:stage["status"]="completed"
        save_run_state(root,state);return root
    def config(self,base:Path)->Path:
        env=base/"boltz2/bin";env.mkdir(parents=True);exe=env/"boltz"
        exe.write_text("""#!/usr/bin/env python3
import json,pathlib,sys
inp=pathlib.Path(sys.argv[2]);out=pathlib.Path(sys.argv[sys.argv.index('--out_dir')+1]);dest=out/'predictions'/inp.stem;dest.mkdir(parents=True,exist_ok=True)
score={'reference':.5,'design_a':.6,'design_b':.5}[inp.stem]
(dest/f'affinity_{inp.stem}.json').write_text(json.dumps({'affinity_probability_binary':score,'affinity_pred_value':1.2}))
(dest/f'{inp.stem}_model_0.cif').write_text('data_test')
""");exe.chmod(0o755);config=base/"config.json";config.write_text(json.dumps({"environments":{"boltz2":str(base/"boltz2")},"boltz2":{"devices":1,"accelerator":"gpu","use_msa_server":True}}));return config
    def test_preparation_writes_reference_and_candidate_yaml(self):
        with tempfile.TemporaryDirectory() as d:
            result=prepare_boltz2_affinity(self.project(Path(d)));self.assertEqual(result["record_ids"],["reference","design_a","design_b"]);text=result["files"]["reference"].read_text();self.assertIn("affinity:",text);self.assertIn("binder: B",text)
    def test_strict_reference_comparison_and_full_run(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);root=self.project(base);code=run_boltz2_affinity(root,config_path=self.config(base));self.assertEqual(code,0)
            with (root/"boltz2_affinity/results/boltz2_affinity.csv").open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual(rows[0]["decision"],"advance_to_structural_confidence");self.assertEqual(rows[1]["decision"],"does_not_advance");self.assertEqual((root/"boltz2_affinity/results/affinity_qualified_candidates.fasta").read_text(),">design_a\nHIKLMN\n")
            plot=root/"boltz2_affinity/results/boltz2_affinity.svg";ET.parse(plot);svg=plot.read_text();self.assertIn('class="tick"',svg);self.assertIn('>0.60</text>',svg);self.assertNotIn('>design_a</text>',svg);self.assertIn('<title>design_a:',svg);self.assertIn('advance if probability &gt; 0.500',svg)
            state=load_run_state(root);self.assertEqual(next(s for s in state["stages"] if s["name"]=="boltz2_affinity")["status"],"completed")
    def test_noncatalytic_mode_uses_shortlist_when_liability_is_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            root=self.project(Path(d),"non-catalytic-reference");state=load_run_state(root);next(s for s in state["stages"] if s["name"]=="catalytic_liability")["status"]="skipped";save_run_state(root,state)
            (root/"catalytic_liability/input").mkdir(parents=True);(root/"catalytic_liability/input/shortlisted_candidates.fasta").write_text(">design_a\nHIKLMN\n")
            result=prepare_boltz2_affinity(root);self.assertEqual(result["mode"],"non-catalytic-reference");self.assertEqual(result["candidate_ids"],["design_a"])
    def test_zero_survivors_automatically_skip_confidence_and_write_report(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);root=self.project(base);config=self.config(base);exe=base/"boltz2/bin/boltz"
            exe.write_text(exe.read_text().replace("'design_a':.6", "'design_a':.4"))
            self.assertEqual(run_boltz2_affinity(root,config_path=config),0)
            state=load_run_state(root)
            self.assertEqual(next(s for s in state["stages"] if s["name"]=="boltz2_affinity")["status"],"completed")
            self.assertEqual(next(s for s in state["stages"] if s["name"]=="structural_confidence")["status"],"skipped")
            self.assertEqual(next(s for s in state["stages"] if s["name"]=="reports")["status"],"completed")
            self.assertTrue((root/"reports/terminal_summary.md").is_file())
if __name__=="__main__":unittest.main()
