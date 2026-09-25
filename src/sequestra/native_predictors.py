"""Standalone native DLKcat and CatPred runners (stdlib only)."""
from __future__ import annotations
import argparse,csv,math,os,shutil,subprocess,sys
from pathlib import Path
from typing import Optional

def _rows(path:Path)->list[dict[str,str]]:
    with path.open(newline="",encoding="utf-8-sig") as handle:return list(csv.DictReader(handle))
def _write(path:Path,fields:list[str],rows:list[dict[str,object]],delimiter:str=",")->None:
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="",encoding="utf-8") as handle:w=csv.DictWriter(handle,fieldnames=fields,delimiter=delimiter);w.writeheader();w.writerows(rows)
def _run(command:list[str],cwd:Path)->None:
    print("NATIVE COMMAND:",subprocess.list2cmdline(command),flush=True)
    completed=subprocess.run(command,cwd=cwd,check=False)
    if completed.returncode:raise RuntimeError(f"native predictor exited with status {completed.returncode}")
def _sequence_map(rows:list[dict[str,str]])->dict[str,list[str]]:
    result:dict[str,list[str]]={}
    for row in rows:
        sequence="".join(row["sequence"].split()).upper()
        result.setdefault(sequence,[]).append(row["candidate_id"])
    return result
def _take_identity(mapping:dict[str,list[str]],sequence:str)->Optional[str]:
    identifiers=mapping.get("".join(sequence.split()).upper(),[])
    return identifiers.pop(0) if identifiers else None
def _expected(mapping:dict[str,list[str]])->set[str]:return {candidate for identifiers in mapping.values() for candidate in identifiers}

def _catpred_identity(row:dict[str,str],mapping:dict[str,list[str]],expected:set[str])->Optional[str]:
    by_sequence=_take_identity(mapping,row.get("sequence",""))
    path_id=Path(row.get("pdbpath","")).stem
    by_path=path_id if path_id in expected else None
    if by_sequence and by_path and by_sequence!=by_path:raise ValueError(f"CatPred identity disagreement: sequence={by_sequence}, pdbpath={by_path}")
    return by_sequence or by_path

def _convert_mm(value_mm:float,unit:str)->float:
    normalized=unit.strip().replace("µ","u").replace("μ","u").lower()
    if normalized=="m":return value_mm/1000.0
    if normalized=="mm":return value_mm
    if normalized=="um":return value_mm*1000.0
    raise ValueError(f"unsupported project Km unit for CatPred conversion: {unit!r}; use M, mM, or uM")

def _normalize_catpred(inputs:list[dict[str,str]],raw:Path,output:Path)->None:
    mapping=_sequence_map(inputs);expected=_expected(mapping);normalized=[]
    units={row.get("km_unit","").strip() or "M" for row in inputs}
    if len(units)!=1:raise ValueError(f"CatPred inputs contain inconsistent Km units: {sorted(units)}")
    unit=next(iter(units))
    for row in _rows(raw):
        cid=_catpred_identity(row,mapping,expected)
        log_value=row.get("log10km_mean","").strip()
        linear=row.get("Prediction_(mM)","").strip() or row.get("prediction_mM","").strip()
        if linear:
            value_mm=float(linear)
            if not log_value and value_mm>0:log_value=str(math.log10(value_mm))
        elif log_value:value_mm=10.0**float(log_value)
        else:raise ValueError("CatPred output lacks both log10km_mean and Prediction_(mM)")
        variance=row.get("log10km_mean_mve_uncal_var","").strip()
        log_sd=str(math.sqrt(max(0.0,float(variance)))) if variance else row.get("SD_total","").strip()
        if cid:normalized.append({"candidate_id":cid,"predicted_km":_convert_mm(value_mm,unit),"predicted_km_unit":unit,"catpred_log10km_mean_mM":log_value,"km_uncertainty":log_sd,"km_uncertainty_scale":"log10_mM_sd"})
    if set(r["candidate_id"] for r in normalized)!=expected:raise ValueError("CatPred output could not be reconciled to every input identity")
    _write(output,["candidate_id","predicted_km","predicted_km_unit","catpred_log10km_mean_mM","km_uncertainty","km_uncertainty_scale"],normalized)

def run_dlkcat(args:argparse.Namespace)->None:
    repository=Path(args.repository).expanduser().resolve();output=Path(args.output_dir).resolve();output.mkdir(parents=True,exist_ok=True)
    script=repository/"DeeplearningApproach/Code/example/prediction_for_input.py"
    data=repository/"DeeplearningApproach/Data";results=repository/"DeeplearningApproach/Results"
    for path in (script,data,results):
        if not path.exists():raise FileNotFoundError(f"required DLKcat path not found: {path}")
    inputs=_rows(Path(args.input_csv));mapping=_sequence_map(inputs);expected=_expected(mapping);native=output/"dlkcat_input.tsv"
    _write(native,["Substrate Name","Substrate SMILES","Protein Sequence"],[{"Substrate Name":r["substrate_name"],"Substrate SMILES":r["substrate_smiles"],"Protein Sequence":r["sequence"]} for r in inputs],"\t")
    stage=output/"native_work/DeeplearningApproach/Code/example";stage.mkdir(parents=True,exist_ok=True)
    staged_root=stage.parents[1]
    for name,target in (("Data",data),("Results",results)):
        link=staged_root/name
        if link.is_symlink() and link.resolve()==target:continue
        if link.exists() or link.is_symlink():raise RuntimeError(f"unsafe existing DLKcat staging path: {link}")
        link.symlink_to(target,target_is_directory=True)
    _run([sys.executable,str(script),str(native)],stage)
    raw=stage/"output.tsv"
    if not raw.is_file():raise FileNotFoundError(f"DLKcat produced no output: {raw}")
    shutil.copy2(raw,output/"dlkcat_native_output.tsv")
    normalized=[]
    for row in _rows_tsv(raw):
        cid=_take_identity(mapping,row.get("Protein Sequence",""))
        if cid:normalized.append({"candidate_id":cid,"predicted_kcat":row.get("Kcat value (1/s)","")})
    if set(r["candidate_id"] for r in normalized)!=expected:raise ValueError("DLKcat output could not be reconciled to every input sequence")
    _write(output/"predictions.csv",["candidate_id","predicted_kcat"],normalized)

def _rows_tsv(path:Path)->list[dict[str,str]]:
    with path.open(newline="",encoding="utf-8-sig") as handle:return list(csv.DictReader(handle,delimiter="\t"))
def _checkpoint(data_root:Path,configured:Optional[str])->Path:
    if configured:
        path=Path(configured).expanduser().resolve()
        if path.is_dir():return path
        raise FileNotFoundError(f"configured CatPred Km checkpoint directory not found: {path}")
    candidates=[p for p in data_root.rglob("*") if p.is_dir() and p.name.lower()=="km" and any(x.is_file() for x in p.rglob("*"))]
    candidates.sort(key=lambda p:("production" not in str(p).lower(),"reproduce_checkpoints" not in str(p).lower(),len(p.parts)))
    if not candidates:raise FileNotFoundError(f"no CatPred Km checkpoint directory found under {data_root}")
    return candidates[0]

def run_catpred(args:argparse.Namespace)->None:
    repository=Path(args.repository).expanduser().resolve();data_root=Path(args.data_root).expanduser().resolve();output=Path(args.output_dir).resolve();output.mkdir(parents=True,exist_ok=True)
    create_records=repository/"scripts/create_pdbrecords.py";predict=repository/"predict.py"
    for path in (repository,data_root,create_records,predict):
        if not path.exists():raise FileNotFoundError(f"required CatPred path not found: {path}")
    inputs=_rows(Path(args.input_csv));mapping=_sequence_map(inputs);expected=_expected(mapping);native=output/"catpred_input.csv";records=output/"catpred_input.json.gz";raw=output/"catpred_native_output.csv"
    _write(native,["Substrate","SMILES","sequence","pdbpath"],[{"Substrate":r["substrate_name"],"SMILES":r["substrate_smiles"],"sequence":r["sequence"],"pdbpath":r["candidate_id"]+".pdb"} for r in inputs])
    checkpoint=_checkpoint(data_root,args.checkpoint_dir)
    if raw.is_file():print(f"Reusing existing CatPred raw predictions: {raw}",flush=True)
    else:
        _run([sys.executable,str(create_records),"--data_file",str(native),"--out_file",str(records)],repository)
        command=[sys.executable,str(predict),"--test_path",str(native),"--preds_path",str(raw),"--checkpoint_dir",str(checkpoint),"--uncertainty_method","mve","--smiles_column","SMILES","--individual_ensemble_predictions","--protein_records_path",str(records)]
        _run(command,repository)
    _normalize_catpred(inputs,raw,output/"predictions.csv")

def main(argv:Optional[list[str]]=None)->int:
    parser=argparse.ArgumentParser();sub=parser.add_subparsers(dest="predictor",required=True)
    for name in ("dlkcat","catpred"):
        p=sub.add_parser(name);p.add_argument("--repository",required=True);p.add_argument("--input-csv",required=True);p.add_argument("--output-dir",required=True)
        if name=="catpred":p.add_argument("--data-root",required=True);p.add_argument("--checkpoint-dir");p.add_argument("--use-gpu",action="store_true")
    args=parser.parse_args(argv);run_dlkcat(args) if args.predictor=="dlkcat" else run_catpred(args);return 0
if __name__=="__main__":raise SystemExit(main())
