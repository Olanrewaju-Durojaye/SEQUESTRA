"""Direct, checkpointed DLKcat/CatPred catalytic-liability integration."""
from __future__ import annotations
import csv, html, json, math, re, shlex, subprocess
from pathlib import Path
from typing import Any
from .decision_rules import CATALYTIC_REFERENCE, classify_catalytic_liability
from .reference import read_fasta
from .run_state import finish_stage, record_partial_stage_success, start_stage

SEQUENCE_COLUMNS=("designed_chain_sequence","designed_sequence","sequence","protein_sequence")
ID_COLUMNS=("candidate_id","design_id","id","name","sequestra_candidate_id")
KCAT_COLUMNS=("kcat","predicted_kcat","kcat_prediction","kcat_s","kcat_s^-1")
KM_COLUMNS=("km","predicted_km","km_prediction","km_m","km_mm")
UNCERTAINTY_COLUMNS=("km_uncertainty","catpred_uncertainty","uncertainty","km_sd_total")

def _load(path: Path)->dict[str,Any]: return json.loads(path.read_text(encoding="utf-8"))
def _first(row:dict[str,str], cols:tuple[str,...])->str|None:
    low={k.lower():v for k,v in row.items()}
    for col in cols:
        value=low.get(col.lower())
        if value is not None and str(value).strip() not in {"","NA","N/A","null","None"}: return str(value).strip()
    return None
def _number(value:str|None)->float|None:
    try: result=float(value) if value is not None else None
    except ValueError: return None
    return result if result is not None and math.isfinite(result) else None
def _read_csv(path:Path)->list[dict[str,str]]:
    with path.open(newline="",encoding="utf-8-sig") as handle: return list(csv.DictReader(handle))

def _sequence_from_cif(path:Path)->str:
    if not path.is_file(): return ""
    text=path.read_text(encoding="utf-8",errors="replace")
    for key in ("_entity_poly.pdbx_seq_one_letter_code_can","_entity_poly.pdbx_seq_one_letter_code"):
        match=re.search(rf"(?ms)^{re.escape(key)}\s*\n;\s*(.*?)\s*;",text)
        if match:
            seq=re.sub(r"[^A-Za-z]","",match.group(1)).upper()
            if seq:return seq
    return ""

def prepare_catalytic_liability(project_dir:Path)->dict[str,Any]:
    root=project_dir.expanduser().resolve(); spec=_load(root/"configuration/project_specification.yaml")
    if not spec["decision_rules"]["catalytic_liability_enabled"]: raise ValueError("catalytic-liability screening is disabled")
    selection=root/"shortlist/selection/selected_candidates.csv"
    if not selection.is_file(): raise FileNotFoundError(f"shortlist selection not found: {selection}")
    rows=_read_csv(selection)
    if not rows: raise ValueError("shortlist contains no candidates")
    ref_records=read_fasta(root/"inputs/reference_complete.fasta")
    if len(ref_records)!=1: raise ValueError("reference FASTA must contain exactly one sequence")
    records=[("reference",ref_records[0].sequence)]; seen={"reference"}
    for row in rows:
        cid=row["sequestra_candidate_id"].strip()
        if not cid or cid in seen: raise ValueError(f"invalid or duplicate candidate identifier: {cid!r}")
        seen.add(cid); raw=next((row.get(k,"") for k in SEQUENCE_COLUMNS if row.get(k,"").strip()),"")
        seq="".join(raw.split()).upper() or _sequence_from_cif(root/"shortlist/selection/selected_binders"/(row.get("file_name") or ""))
        if not seq or any(a not in "ACDEFGHIKLMNPQRSTVWY" for a in seq): raise ValueError(f"invalid or missing sequence for {cid}")
        records.append((cid,seq))
    work=root/"catalytic_liability/input"; work.mkdir(parents=True,exist_ok=True)
    candidates=work/"shortlisted_candidates.fasta"; candidates.write_text("".join(f">{i}\n{s}\n" for i,s in records[1:]),encoding="utf-8")
    all_fasta=work/"reference_and_candidates.fasta"; all_fasta.write_text("".join(f">{i}\n{s}\n" for i,s in records),encoding="utf-8")
    input_csv=work/"predictor_input.csv"
    km_unit=spec["decision_rules"].get("km_unit") or "M"
    with input_csv.open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=["candidate_id","sequence","substrate_name","substrate_smiles","km_unit"]); writer.writeheader()
        for cid,seq in records: writer.writerow({"candidate_id":cid,"sequence":seq,"substrate_name":spec["ligand"]["name"],"substrate_smiles":spec["ligand"]["smiles"],"km_unit":km_unit})
    manifest={"schema_version":2,"reference_mode":spec["reference"]["mode"],"candidate_count":len(records)-1,"candidate_ids":[x[0] for x in records[1:]],"record_ids":[x[0] for x in records],"substrate_name":spec["ligand"]["name"],"substrate_smiles":spec["ligand"]["smiles"],"kcat_unit":spec["decision_rules"].get("kcat_unit") or "s^-1","km_unit":km_unit,"prediction_failure_establishes_non_catalysis":False}
    manifest_path=work/"input_manifest.json"; manifest_path.write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
    return {**manifest,"manifest_path":manifest_path,"candidate_fasta_path":candidates,"all_fasta_path":all_fasta,"predictor_input_path":input_csv}

def _legacy(root:Path)->tuple[Path,dict[str,dict[str,str]]]|None:
    old=root/"catalytic_liability/enztra_job"
    if not old.is_dir(): return None
    fallback=None
    for path in sorted(old.rglob("*.csv"),key=lambda p:(p.name!="kinetics.csv",len(p.parts))):
        rows=_read_csv(path); fields={k.lower() for k in (rows[0] if rows else {})}
        if any(k.lower() in fields for k in KCAT_COLUMNS) and any(k.lower() in fields for k in KM_COLUMNS):
            mapped={i:r for r in rows if (i:=_first(r,ID_COLUMNS))}
            if any((_first(row,("role","record_type","type")) or "").lower()=="reference" for row in rows):
                return path,mapped
            fallback=(path,mapped)
    return fallback

def _adapter(config:dict[str,Any],name:str,root:Path,prepared:dict[str,Any])->tuple[list[str],Path]:
    cfg=config.get("predictors",{}).get(name)
    if not isinstance(cfg,dict): raise KeyError(f"configuration lacks predictors.{name}; see sequestra.config.json")
    env=Path(config["environments"][name]).expanduser().resolve(); repo=Path(cfg["repository"]).expanduser().resolve(); outdir=root/f"catalytic_liability/{name}"; outdir.mkdir(parents=True,exist_ok=True)
    values={"python":str(env/"bin/python"),"environment":str(env),"repository":str(repo),"project":str(root),"input_csv":str(prepared["predictor_input_path"]),"input_fasta":str(prepared["all_fasta_path"]),"output_dir":str(outdir),"substrate_name":prepared["substrate_name"],"substrate_smiles":prepared["substrate_smiles"]}
    if cfg.get("adapter")=="native":
        if not Path(values["python"]).is_file():raise FileNotFoundError(f"{name} environment Python not found: {values['python']}")
        runner=Path(__file__).with_name("native_predictors.py").resolve()
        command=[values["python"],str(runner),name,"--repository",str(repo),"--input-csv",values["input_csv"],"--output-dir",str(outdir)]
        if name=="dlkcat":
            required=[repo/"DeeplearningApproach/Code/example/prediction_for_input.py",repo/"DeeplearningApproach/Data",repo/"DeeplearningApproach/Results"]
            missing=[str(path) for path in required if not path.exists()]
            if missing:raise FileNotFoundError("missing native DLKcat resources: "+", ".join(missing))
        else:
            data_root=Path(cfg["data_root"]).expanduser().resolve();required=[repo/"predict.py",repo/"scripts/create_pdbrecords.py",data_root]
            missing=[str(path) for path in required if not path.exists()]
            if missing:raise FileNotFoundError("missing native CatPred resources: "+", ".join(missing))
            checkpoint=Path(cfg["checkpoint_dir"]).expanduser().resolve() if cfg.get("checkpoint_dir") else next((path for path in sorted(data_root.rglob("km"),key=lambda p:("production" not in str(p).lower(),len(p.parts))) if path.is_dir() and any(item.is_file() for item in path.rglob("*"))),None)
            if checkpoint is None or not checkpoint.is_dir():raise FileNotFoundError(f"no CatPred Km checkpoint directory found under {data_root}")
            command += ["--data-root",str(data_root),"--checkpoint-dir",str(checkpoint)]
            if cfg.get("use_gpu",True):command.append("--use-gpu")
        return command,outdir/"predictions.csv"
    template=cfg.get("command")
    if not isinstance(template,list) or not template: raise ValueError(f"predictors.{name}.command must be a non-empty array")
    return [str(x).format(**values) for x in template],Path(str(cfg["output_csv"]).format(**values)).expanduser()

def build_predictor_commands(project_dir:Path,*,sequestra_config_path:Path)->dict[str,tuple[list[str],Path]]:
    root=project_dir.expanduser().resolve(); prepared=prepare_catalytic_liability(root); config=_load(sequestra_config_path.expanduser().resolve())
    return {name:_adapter(config,name,root,prepared) for name in ("dlkcat","catpred")}

def _checkpoint(root:Path,name:str)->Path:return root/f"catalytic_liability/checkpoints/{name}.json"
def _mark(root:Path,name:str,status:str,**extra:Any)->None:
    path=_checkpoint(root,name); path.parent.mkdir(parents=True,exist_ok=True); path.write_text(json.dumps({"schema_version":1,"status":status,**extra},indent=2)+"\n",encoding="utf-8")
def _complete(root:Path,name:str,output:Path)->bool:
    try:
        if not output.is_file() or _load(_checkpoint(root,name)).get("status")!="completed":return False
        if name=="catpred":
            rows=_read_csv(output)
            return bool(rows) and all(_first(row,KM_COLUMNS) is not None for row in rows)
        return True
    except (OSError,json.JSONDecodeError):return False
def _prediction_map(path:Path,columns:tuple[str,...])->dict[str,dict[str,str]]:
    if not path.is_file():raise FileNotFoundError(f"prediction table not found: {path}")
    return {i:r for r in _read_csv(path) if (i:=_first(r,ID_COLUMNS)) and _first(r,columns) is not None}

def _classify_results(root:Path,dlkcat_table:Path|None=None,catpred_table:Path|None=None)->dict[str,Any]:
    prepared=prepare_catalytic_liability(root); old=_legacy(root) if dlkcat_table is None or catpred_table is None else None
    if old:
        sources=[old[0]]; krows=mrows=old[1]
        if "reference" not in krows:
            role_reference=next((key for key,row in krows.items() if (_first(row,("role","record_type","type")) or "").lower()=="reference"),None)
            unmatched=[key for key in krows if key not in set(prepared["candidate_ids"])]
            reference_key=role_reference or (unmatched[0] if len(unmatched)==1 else None)
            if reference_key:
                krows["reference"]=mrows["reference"]=krows[reference_key]
    else:
        if dlkcat_table is None or catpred_table is None:raise FileNotFoundError("both direct predictor results are required")
        sources=[dlkcat_table,catpred_table]; krows=_prediction_map(dlkcat_table,KCAT_COLUMNS); mrows=_prediction_map(catpred_table,KM_COLUMNS)
    expected=set(prepared["record_ids"]); mk=sorted(expected-set(krows)); mm=sorted(expected-set(mrows))
    if mk or mm:raise ValueError(f"predictor identity validation failed; missing DLKcat={mk}, CatPred={mm}")
    rk=_number(_first(krows["reference"],KCAT_COLUMNS)); rm=_number(_first(mrows["reference"],KM_COLUMNS)); informational=prepared["reference_mode"]!=CATALYTIC_REFERENCE; output=[]
    for cid in prepared["candidate_ids"]:
        k=_number(_first(krows[cid],KCAT_COLUMNS)); m=_number(_first(mrows[cid],KM_COLUMNS)); result=classify_catalytic_liability(candidate_kcat=k,candidate_km=m,reference_kcat=rk,reference_km=rm,within_supported_interpretation=not informational)
        output.append({"candidate_id":cid,"predicted_kcat":k,"kcat_unit":prepared["kcat_unit"],"predicted_km":m,"km_unit":prepared["km_unit"],"km_uncertainty_log10_sd":_number(_first(mrows[cid],UNCERTAINTY_COLUMNS)),"reference_kcat":rk,"reference_km":rm,"classification":result["classification"],"detail":result["detail"],"decision":"informational_only_no_exclusion" if informational else result["decision"]})
    dest=root/"catalytic_liability/results"; dest.mkdir(parents=True,exist_ok=True); result_path=dest/"catalytic_liability.csv"
    with result_path.open("w",newline="",encoding="utf-8") as handle:w=csv.DictWriter(handle,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    survivor_decisions={"advance","advance_with_partial_liability_flag","informational_only_no_exclusion"}; survivors=[r for r in output if r["decision"] in survivor_decisions]; review=[r for r in output if r["decision"] not in survivor_decisions|{"exclude"}]
    seqs={r.header.split()[0]:r.sequence for r in read_fasta(root/"catalytic_liability/input/shortlisted_candidates.fasta")}; survivor_path=dest/"automatic_survivors.fasta"; survivor_path.write_text("".join(f">{r['candidate_id']}\n{seqs[r['candidate_id']]}\n" for r in survivors),encoding="utf-8")
    plot_path=dest/"catalytic_liability.svg"; _write_plot(plot_path,output,rk,rm)
    summary={"schema_version":2,"source_tables":[str(p.relative_to(root)) for p in sources],"reference_kcat":rk,"reference_km":rm,"kcat_unit":prepared["kcat_unit"],"km_unit":prepared["km_unit"],"candidate_count":len(output),"excluded_count":sum(r["decision"]=="exclude" for r in output),"automatic_survivor_count":len(survivors),"manual_review_count":len(review),"automatic_survivor_ids":[r["candidate_id"] for r in survivors],"manual_review_ids":[r["candidate_id"] for r in review],"interpretation_warning":"Predictions are hypotheses; missing predictions are indeterminate, never evidence of non-catalysis."}; summary_path=dest/"summary.json";summary_path.write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    return {**summary,"result_path":result_path,"summary_path":summary_path,"survivor_path":survivor_path,"plot_path":plot_path}

def _write_plot(path:Path,rows:list[dict[str,Any]],rk:float|None,rm:float|None)->None:
    valid=[r for r in rows if r["predicted_kcat"] and r["predicted_km"] and r["predicted_kcat"]>0 and r["predicted_km"]>0]; pts=[(math.log10(r["predicted_km"]),math.log10(r["predicted_kcat"]),r) for r in valid]
    if rk and rm and rk>0 and rm>0:pts.append((math.log10(rm),math.log10(rk),{"candidate_id":"reference","decision":"reference"}))
    xs=[p[0] for p in pts] or [-1,1];ys=[p[1] for p in pts] or [-1,1];xmin,xmax=min(xs)-.2,max(xs)+.2;ymin,ymax=min(ys)-.2,max(ys)+.2
    if xmin==xmax:xmin-=.5;xmax+=.5
    if ymin==ymax:ymin-=.5;ymax+=.5
    sx=lambda x:90+(x-xmin)/(xmax-xmin)*640;sy=lambda y:450-(y-ymin)/(ymax-ymin)*400
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="760" height="540" viewBox="0 0 760 540">','<rect width="100%" height="100%" fill="white"/>','<style>text{font-family:sans-serif;font-size:12px}.title{font-size:18px;font-weight:bold}.tick{font-size:11px;fill:#444}.grid{stroke:#e0e0e0;stroke-width:1}</style>','<text class="title" x="90" y="28">SEQUESTRA exclusion-oriented catalytic-liability map</text>']
    if rk and rm and rk>0 and rm>0:
        rx,ry=sx(math.log10(rm)),sy(math.log10(rk));parts += [f'<rect x="90" y="50" width="{max(0,rx-90):.1f}" height="{max(0,ry-50):.1f}" fill="#ffebee"/>']
    for index in range(5):
        value=xmin+(xmax-xmin)*index/4; x=sx(value); parts += [f'<line class="grid" x1="{x:.1f}" y1="50" x2="{x:.1f}" y2="450"/>',f'<line x1="{x:.1f}" y1="450" x2="{x:.1f}" y2="456" stroke="black"/>',f'<text class="tick" x="{x:.1f}" y="472" text-anchor="middle">{value:.2f}</text>']
    for index in range(5):
        value=ymin+(ymax-ymin)*index/4; y=sy(value); parts += [f'<line class="grid" x1="90" y1="{y:.1f}" x2="730" y2="{y:.1f}"/>',f'<line x1="84" y1="{y:.1f}" x2="90" y2="{y:.1f}" stroke="black"/>',f'<text class="tick" x="78" y="{y+4:.1f}" text-anchor="end">{value:.2f}</text>']
    parts += ['<line x1="90" y1="450" x2="730" y2="450" stroke="black"/><line x1="90" y1="50" x2="90" y2="450" stroke="black"/>']
    if rk and rm and rk>0 and rm>0:
        parts += [f'<line x1="{rx:.1f}" y1="50" x2="{rx:.1f}" y2="450" stroke="#777" stroke-dasharray="5 4"/>',f'<line x1="90" y1="{ry:.1f}" x2="730" y2="{ry:.1f}" stroke="#777" stroke-dasharray="5 4"/>','<text x="98" y="68" fill="#b71c1c">EXCLUDE: higher kcat + lower Km</text>']
    colors={"exclude":"#c62828","advance":"#2e7d32","advance_with_partial_liability_flag":"#ef6c00","reference":"#1565c0","informational_only_no_exclusion":"#6a1b9a"}
    for x,y,r in pts:
        label=html.escape(str(r["candidate_id"])); decision=html.escape(str(r["decision"])); parts += [f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="6" fill="{colors.get(r["decision"],"#616161")}"><title>{label}: log10 Km={x:.3f}, log10 kcat={y:.3f}; {decision}</title></circle>']
    parts += ['<text x="310" y="505">log10 predicted Km</text>','<text transform="translate(22 315) rotate(-90)">log10 predicted kcat</text>','<text x="90" y="530">Red: excluded | Orange: partial flag | Green: advance | Blue: reference</text>','</svg>'];path.write_text("\n".join(parts)+"\n",encoding="utf-8")

def run_catalytic_liability(project_dir:Path,*,sequestra_config_path:Path,prepare_only:bool=False)->int:
    root=project_dir.expanduser().resolve();prepared=prepare_catalytic_liability(root);old=_legacy(root);commands={} if old else build_predictor_commands(root,sequestra_config_path=sequestra_config_path);preview=[["import-existing",str(old[0])]] if old else [commands[n][0] for n in ("dlkcat","catpred")];start_stage(root,"catalytic_liability",command=[x for cmd in preview for x in cmd]);log_path=root/"logs/catalytic_liability.log";log_path.parent.mkdir(parents=True,exist_ok=True)
    try:
        if prepare_only:record_partial_stage_success(root,"catalytic_liability",message="direct predictor inputs prepared; inference not executed",outputs=[str(prepared["manifest_path"].relative_to(root))]);print("Preparation completed. Stage remains resumable.");return 0
        if old:
            _mark(root,"dlkcat","completed",source="legacy_combined",output=str(old[0].relative_to(root)));_mark(root,"catpred","completed",source="legacy_combined",output=str(old[0].relative_to(root)));result=_classify_results(root)
        else:
            with log_path.open("a",encoding="utf-8") as log:
                for name in ("dlkcat","catpred"):
                    command,output=commands[name]
                    if _complete(root,name,output):log.write(f"\nRESUME: reusing completed {name}: {output}\n");continue
                    _mark(root,name,"running",command=command,output=str(output));log.write("\nCOMMAND: "+shlex.join(command)+"\n");log.flush();proc=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                    assert proc.stdout is not None
                    for line in proc.stdout:print(line,end="");log.write(line);log.flush()
                    if proc.wait():raise RuntimeError(f"{name} predictor failed")
                    if not output.is_file():raise FileNotFoundError(f"{name} configured output is absent: {output}")
                    _mark(root,name,"completed",command=command,output=str(output))
            result=_classify_results(root,commands["dlkcat"][1],commands["catpred"][1])
        _mark(root,"combination","completed",output=str(result["result_path"].relative_to(root)));outputs=[str(result[k].relative_to(root)) for k in ("result_path","summary_path","survivor_path","plot_path")];finish_stage(root,"catalytic_liability",succeeded=True,message=f"classified {result['candidate_count']} candidates; {result['excluded_count']} excluded",outputs=outputs);print(f"Catalytic-liability completed: {result['excluded_count']} excluded; {result['automatic_survivor_count']} advance; {result['manual_review_count']} review.");return 0
    except Exception as error:finish_stage(root,"catalytic_liability",succeeded=False,message=str(error));raise
