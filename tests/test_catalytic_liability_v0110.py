from __future__ import annotations
import csv,json,tempfile,unittest
from pathlib import Path
from sequestra.catalytic_liability import _classify_results,_complete,_mark,build_predictor_commands,prepare_catalytic_liability,run_catalytic_liability

class CatalyticLiabilityTests(unittest.TestCase):
    def test_completed_catpred_checkpoint_with_blank_values_is_not_reused(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);output=root/"catalytic_liability/catpred/predictions.csv";output.parent.mkdir(parents=True)
            output.write_text("candidate_id,predicted_km\nreference,\n")
            _mark(root,"catpred","completed",output=str(output))
            self.assertFalse(_complete(root,"catpred",output))
    def project(self,base:Path)->Path:
        root=base/"project"
        for p in (root/"configuration",root/"inputs",root/"shortlist/selection",root/"catalytic_liability/enztra_job"):p.mkdir(parents=True,exist_ok=True)
        spec={"reference":{"mode":"catalytic-reference"},"ligand":{"name":"PMM","smiles":"CCO"},"decision_rules":{"catalytic_liability_enabled":True,"kcat_unit":"s^-1","km_unit":"M"}}
        (root/"configuration/project_specification.yaml").write_text(json.dumps(spec));(root/"inputs/reference_complete.fasta").write_text(">ref\nACDEFGHIK\n")
        rows=[{"sequestra_candidate_id":"design_a","sequence":"ACDEFGHIK","file_name":"a.cif"},{"sequestra_candidate_id":"design_b","sequence":"LMNPQRSTV","file_name":"b.cif"}]
        with (root/"shortlist/selection/selected_candidates.csv").open("w",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        return root
    def test_prepares_reference_and_candidates(self):
        with tempfile.TemporaryDirectory() as d:
            result=prepare_catalytic_liability(self.project(Path(d)));self.assertEqual(result["record_ids"],["reference","design_a","design_b"]);self.assertTrue(result["predictor_input_path"].is_file())
    def test_direct_commands_are_independent(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);root=self.project(base);env=base/"env/bin";env.mkdir(parents=True);(env/"python").write_text("");repo=base/"repo";repo.mkdir()
            config=base/"config.json";config.write_text(json.dumps({"environments":{"dlkcat":str(base/"env"),"catpred":str(base/"env")},"predictors":{n:{"repository":str(repo),"command":["{python}","{repository}/predict.py","{input_csv}","{output_dir}"],"output_csv":"{output_dir}/predictions.csv"} for n in ("dlkcat","catpred")}}))
            commands=build_predictor_commands(root,sequestra_config_path=config);self.assertIn("dlkcat",commands);self.assertIn("catpred",commands);self.assertNotEqual(commands["dlkcat"][1],commands["catpred"][1])
    def test_native_commands_validate_real_layouts_before_start(self):
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);root=self.project(base);env=base/"env/bin";env.mkdir(parents=True);(env/"python").write_text("")
            dlk=base/"DLKcat";(dlk/"DeeplearningApproach/Code/example").mkdir(parents=True);(dlk/"DeeplearningApproach/Code/example/prediction_for_input.py").write_text("");(dlk/"DeeplearningApproach/Data").mkdir();(dlk/"DeeplearningApproach/Results").mkdir()
            cat=base/"CatPred";(cat/"scripts").mkdir(parents=True);(cat/"predict.py").write_text("");(cat/"scripts/create_pdbrecords.py").write_text("");data=base/"data/pretrained/production/km";data.mkdir(parents=True);(data/"model.pt").write_text("")
            config=base/"config.json";config.write_text(json.dumps({"environments":{"dlkcat":str(base/"env"),"catpred":str(base/"env")},"predictors":{"dlkcat":{"adapter":"native","repository":str(dlk)},"catpred":{"adapter":"native","repository":str(cat),"data_root":str(base/"data"),"use_gpu":True}}}))
            commands=build_predictor_commands(root,sequestra_config_path=config);self.assertIn("native_predictors.py",commands["dlkcat"][0][1]);self.assertIn("--checkpoint-dir",commands["catpred"][0]);self.assertEqual(commands["catpred"][1],root/"catalytic_liability/catpred/predictions.csv")
    def test_v080_predictions_resume_and_partial_signals_advance(self):
        with tempfile.TemporaryDirectory() as d:
            root=self.project(Path(d));rows=[{"candidate_id":"original_DHPS","kcat":"10","km":"0.1"},{"candidate_id":"design_a","kcat":"11","km":"0.09"},{"candidate_id":"design_b","kcat":"11","km":"0.11"}]
            with (root/"catalytic_liability/enztra_job/selection.csv").open("w",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            result=_classify_results(root);self.assertEqual(result["excluded_count"],1);self.assertEqual(result["automatic_survivor_ids"],["design_b"]);self.assertTrue(result["plot_path"].is_file())
            svg=result["plot_path"].read_text();self.assertIn('class="tick"',svg);self.assertNotIn('>design_a</text>',svg);self.assertIn('<title>design_a:',svg)
    def test_actual_enztra_kinetics_columns_and_role_are_imported(self):
        with tempfile.TemporaryDirectory() as d:
            root=self.project(Path(d)); rows=[{"design_id":"2VEG_1|Chains","role":"reference","kcat_s":"0.5249","km_mm":"0.0016328","km_sd_total":"0.4159"},{"design_id":"design_a","role":"design","kcat_s":"668.6918","km_mm":"0.0234566","km_sd_total":"0.9275"},{"design_id":"design_b","role":"design","kcat_s":"56.2709","km_mm":"0.0128782","km_sd_total":"0.8100"}]
            with (root/"catalytic_liability/enztra_job/kinetics.csv").open("w",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            result=_classify_results(root);self.assertEqual(result["excluded_count"],0);self.assertEqual(result["automatic_survivor_ids"],["design_a","design_b"]);self.assertEqual(result["source_tables"],["catalytic_liability/enztra_job/kinetics.csv"])
if __name__=="__main__":unittest.main()
