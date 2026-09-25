from __future__ import annotations
import csv,tempfile,unittest
from argparse import Namespace
from pathlib import Path
from sequestra.native_predictors import _normalize_catpred,run_catpred,run_dlkcat

class NativePredictorTests(unittest.TestCase):
    def input_csv(self,root:Path)->Path:
        path=root/"input.csv";rows=[{"candidate_id":"reference","sequence":"ACDE","substrate_name":"PMM","substrate_smiles":"CCO"},{"candidate_id":"design_a","sequence":"FGHI","substrate_name":"PMM","substrate_smiles":"CCO"}]
        with path.open("w",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
        return path
    def test_dlkcat_native_adapter_normalizes_by_sequence(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);repo=root/"DLKcat";example=repo/"DeeplearningApproach/Code/example";example.mkdir(parents=True);(repo/"DeeplearningApproach/Data").mkdir();(repo/"DeeplearningApproach/Results").mkdir()
            (example/"prediction_for_input.py").write_text("""import csv,sys\nrows=list(csv.DictReader(open(sys.argv[1]),delimiter='\\t'))\nwith open('output.tsv','w',newline='') as h:\n w=csv.DictWriter(h,fieldnames=['Substrate Name','Substrate SMILES','Protein Sequence','Kcat value (1/s)'],delimiter='\\t');w.writeheader()\n for i,r in enumerate(rows):w.writerow({**r,'Kcat value (1/s)':str(i+1)})\n""")
            out=root/"out";run_dlkcat(Namespace(repository=str(repo),input_csv=str(self.input_csv(root)),output_dir=str(out)))
            with (out/"predictions.csv").open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual([r["candidate_id"] for r in rows],["reference","design_a"]);self.assertEqual(rows[1]["predicted_kcat"],"2")
    def test_catpred_native_adapter_runs_records_and_prediction(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);repo=root/"CatPred";(repo/"scripts").mkdir(parents=True);data=root/"data/pretrained/production/km";data.mkdir(parents=True);(data/"model.pt").write_text("model")
            (repo/"scripts/create_pdbrecords.py").write_text("""import argparse,pathlib\np=argparse.ArgumentParser();p.add_argument('--data_file');p.add_argument('--out_file');a=p.parse_args();pathlib.Path(a.out_file).write_bytes(b'records')\n""")
            (repo/"predict.py").write_text("""import argparse,csv\np=argparse.ArgumentParser();p.add_argument('--test_path');p.add_argument('--preds_path');p.add_argument('--checkpoint_dir');p.add_argument('--uncertainty_method');p.add_argument('--smiles_column');p.add_argument('--individual_ensemble_predictions',action='store_true');p.add_argument('--protein_records_path');a=p.parse_args();rows=list(csv.DictReader(open(a.test_path)))\nwith open(a.preds_path,'w',newline='') as h:\n w=csv.DictWriter(h,fieldnames=list(rows[0])+['Prediction_(mM)','SD_total']);w.writeheader()\n for i,r in enumerate(rows):w.writerow({**r,'Prediction_(mM)':str(.1+i),'SD_total':str(.01+i)})\n""")
            out=root/"out";run_catpred(Namespace(repository=str(repo),data_root=str(root/"data"),input_csv=str(self.input_csv(root)),output_dir=str(out),checkpoint_dir=None,use_gpu=True))
            with (out/"predictions.csv").open() as handle:rows=list(csv.DictReader(handle))
            self.assertEqual([r["candidate_id"] for r in rows],["reference","design_a"]);self.assertEqual(float(rows[0]["predicted_km"]),0.0001);self.assertEqual(rows[0]["predicted_km_unit"],"M")
    def test_catpred_native_log_schema_is_converted_and_identified(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);inputs=[{"candidate_id":"reference","sequence":"ACDE","km_unit":"M"},{"candidate_id":"design_a","sequence":"FGHI","km_unit":"M"}];raw=root/"raw.csv";output=root/"predictions.csv"
            rows=[{"sequence":"ACDE","pdbpath":"reference.pdb","log10km_mean":"-3","log10km_mean_mve_uncal_var":"0.04"},{"sequence":"FGHI","pdbpath":"design_a.pdb","log10km_mean":"-2","log10km_mean_mve_uncal_var":"0.09"}]
            with raw.open("w",newline="") as h:w=csv.DictWriter(h,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
            _normalize_catpred(inputs,raw,output)
            with output.open() as h:normalized=list(csv.DictReader(h))
            self.assertEqual([r["candidate_id"] for r in normalized],["reference","design_a"]);self.assertEqual(float(normalized[0]["predicted_km"]),1e-6);self.assertEqual(float(normalized[0]["km_uncertainty"]),0.2)
if __name__=="__main__":unittest.main()
