"""Direct, reference-relative Boltz-2 affinity stage."""
from __future__ import annotations
import csv,html,json,math,shlex,subprocess
from pathlib import Path
from typing import Any
from .reference import read_fasta
from .run_state import finish_stage,load_run_state,record_partial_stage_success,start_stage

def _load(path:Path)->dict[str,Any]:return json.loads(path.read_text(encoding="utf-8"))
def _write_json(path:Path,data:dict[str,Any])->None:path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(data,indent=2)+"\n",encoding="utf-8")
def _candidate_fasta(root:Path)->Path:
    state=load_run_state(root);liability=next(row for row in state["stages"] if row["name"]=="catalytic_liability")
    survivor=root/"catalytic_liability/results/automatic_survivors.fasta"
    shortlist=root/"catalytic_liability/input/shortlisted_candidates.fasta"
    if liability["status"]=="completed" and survivor.is_file():return survivor
    if liability["status"]=="skipped" and shortlist.is_file():return shortlist
    fallback=root/"shortlist/selection/selected_candidates.csv"
    if liability["status"]=="skipped" and fallback.is_file():
        rows=list(csv.DictReader(fallback.open(newline="",encoding="utf-8-sig")));generated=root/"boltz2_affinity/input/shortlisted_candidates.fasta";generated.parent.mkdir(parents=True,exist_ok=True)
        generated.write_text("".join(f">{r['sequestra_candidate_id']}\n{''.join((r.get('sequence') or r.get('designed_sequence') or r.get('designed_chain_sequence') or '').split())}\n" for r in rows),encoding="utf-8");return generated
    raise ValueError("no affinity-eligible candidate FASTA is available from the preceding stage")
def _yaml(sequence:str,smiles:str,ccd:str)->str:
    ligand=(f"      smiles: {json.dumps(smiles)}" if smiles else f"      ccd: {json.dumps(ccd)}")
    return f"version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: {json.dumps(sequence)}\n  - ligand:\n      id: B\n{ligand}\nproperties:\n  - affinity:\n      binder: B\n"

def prepare_boltz2_affinity(project_dir:Path)->dict[str,Any]:
    root=project_dir.expanduser().resolve();spec=_load(root/"configuration/project_specification.yaml");records=[("reference",read_fasta(root/"inputs/reference_complete.fasta")[0].sequence)]+[(r.header.split()[0],r.sequence) for r in read_fasta(_candidate_fasta(root))]
    if len(records)<2:raise ValueError("Boltz-2 affinity stage has no candidates")
    ids=[r[0] for r in records]
    if len(ids)!=len(set(ids)):raise ValueError("duplicate Boltz-2 record identifier")
    smiles=str(spec["ligand"].get("smiles") or "").strip();ccd=str(spec["ligand"].get("ccd") or "").strip()
    if not smiles and not ccd:raise ValueError("ligand requires SMILES or CCD for Boltz-2")
    inputs=root/"boltz2_affinity/input";inputs.mkdir(parents=True,exist_ok=True)
    files={}
    for identifier,sequence in records:
        path=inputs/f"{identifier}.yaml";path.write_text(_yaml(sequence,smiles,ccd),encoding="utf-8");files[identifier]=path
    manifest={"schema_version":1,"mode":spec["reference"]["mode"],"metric":"affinity_probability_binary","direction":"higher_is_better","strict_comparison":True,"comparison_margin":float(spec["decision_rules"].get("affinity_comparison_margin",0.0)),"record_ids":ids,"candidate_ids":ids[1:],"ligand_name":spec["ligand"]["name"],"input_files":{k:str(v.relative_to(root)) for k,v in files.items()}}
    path=inputs/"input_manifest.json";_write_json(path,manifest);return {**manifest,"files":files,"manifest_path":path}
def build_boltz2_commands(project_dir:Path,*,config_path:Path)->dict[str,list[str]]:
    root=project_dir.expanduser().resolve();prepared=prepare_boltz2_affinity(root);config=_load(config_path.expanduser().resolve());exe=Path(config["environments"]["boltz2"]).expanduser().resolve()/"bin/boltz"
    if not exe.is_file():raise FileNotFoundError(f"Boltz executable not found: {exe}")
    settings=config.get("boltz2",{});commands={}
    for identifier,path in prepared["files"].items():
        out=root/f"boltz2_affinity/raw/{identifier}";command=[str(exe),"predict",str(path),"--out_dir",str(out),"--devices",str(int(settings.get("devices",1))),"--accelerator",str(settings.get("accelerator","gpu"))]
        if settings.get("use_msa_server",True):command.append("--use_msa_server")
        if settings.get("use_potentials",False):command.append("--use_potentials")
        commands[identifier]=command
    return commands
def _checkpoint(root:Path,identifier:str)->Path:return root/f"boltz2_affinity/checkpoints/{identifier}.json"
def _mark(root:Path,identifier:str,status:str,**extra:Any)->None:_write_json(_checkpoint(root,identifier),{"schema_version":1,"status":status,**extra})
def _find_affinity(output:Path,identifier:str)->Path:
    exact=list(output.rglob(f"affinity_{identifier}.json"));candidates=exact or list(output.rglob("affinity_*.json"))
    if len(candidates)!=1:raise FileNotFoundError(f"expected exactly one Boltz-2 affinity JSON for {identifier}; found {len(candidates)}")
    return candidates[0]
def _completed(root:Path,identifier:str,output:Path)->bool:
    try:return _load(_checkpoint(root,identifier)).get("status")=="completed" and _find_affinity(output,identifier).is_file()
    except (OSError,ValueError,FileNotFoundError,json.JSONDecodeError):return False
def _write_plot(path:Path,rows:list[dict[str,Any]],reference:float,margin:float)->None:
    left,right,top,bottom=90,750,55,450
    count=len(rows);positions=[0]+list(range(1,count+1));xmax=max(1,count)
    sx=lambda x:left+x/xmax*(right-left)
    sy=lambda y:bottom-y*(bottom-top)
    threshold=reference+margin
    parts=['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="550" viewBox="0 0 800 550">','<rect width="100%" height="100%" fill="white"/>','<style>text{font-family:sans-serif;font-size:12px}.title{font-size:18px;font-weight:bold}.tick{font-size:11px;fill:#444}.grid{stroke:#e0e0e0;stroke-width:1}.legend{font-size:12px}</style>','<text class="title" x="90" y="29">SEQUESTRA reference-relative Boltz-2 affinity map</text>']
    if threshold<1:
        shade_top=top if threshold<=0 else sy(threshold)
        parts.append(f'<rect x="{left}" y="{top}" width="{right-left}" height="{max(0,shade_top-top):.1f}" fill="#e8f5e9"/>')
    for index in range(6):
        value=index/5;y=sy(value);parts += [f'<line class="grid" x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}"/>',f'<line x1="{left-6}" y1="{y:.1f}" x2="{left}" y2="{y:.1f}" stroke="black"/>',f'<text class="tick" x="{left-12}" y="{y+4:.1f}" text-anchor="end">{value:.2f}</text>']
    tick_step=max(1,math.ceil(count/10))
    tick_positions=[0]+[i for i in range(1,count+1) if i==1 or i==count or i%tick_step==0]
    for position in sorted(set(tick_positions)):
        x=sx(position);label="R" if position==0 else str(position);parts += [f'<line x1="{x:.1f}" y1="{bottom}" x2="{x:.1f}" y2="{bottom+6}" stroke="black"/>',f'<text class="tick" x="{x:.1f}" y="{bottom+22}" text-anchor="middle">{label}</text>']
    parts += [f'<line x1="{left}" y1="{bottom}" x2="{right}" y2="{bottom}" stroke="black"/>',f'<line x1="{left}" y1="{top}" x2="{left}" y2="{bottom}" stroke="black"/>']
    if 0<=threshold<=1:
        y=sy(threshold);parts += [f'<line x1="{left}" y1="{y:.1f}" x2="{right}" y2="{y:.1f}" stroke="#2e7d32" stroke-width="2" stroke-dasharray="7 5"/>',f'<text x="{right-4}" y="{max(top+14,y-7):.1f}" text-anchor="end" fill="#1b5e20">advance if probability &gt; {threshold:.3f}</text>']
    else:
        parts.append(f'<text x="{right-4}" y="{top+15}" text-anchor="end" fill="#1b5e20">threshold {threshold:.3f} is above plotted probability range</text>')
    points=[(0,"reference",reference,"reference")]+[(i,r["candidate_id"],float(r["affinity_probability_binary"]),r["decision"]) for i,r in enumerate(rows,1)]
    colors={"reference":"#1565c0","advance_to_structural_confidence":"#2e7d32","does_not_advance":"#757575"}
    for position,identifier,value,decision in points:
        label=html.escape(str(identifier));detail=html.escape(str(decision));parts.append(f'<circle cx="{sx(position):.1f}" cy="{sy(value):.1f}" r="6" fill="{colors[decision]}"><title>{label}: affinity_probability_binary={value:.6f}; {detail}</title></circle>')
    parts += [f'<text x="{(left+right)/2:.1f}" y="{bottom+54}" text-anchor="middle">Reference (R) and candidates (input order)</text>','<text transform="translate(24 335) rotate(-90)">Boltz-2 affinity_probability_binary</text>','<circle cx="103" cy="526" r="5" fill="#1565c0"/><text class="legend" x="114" y="530">Reference</text>','<circle cx="205" cy="526" r="5" fill="#2e7d32"/><text class="legend" x="216" y="530">Advance</text>','<circle cx="300" cy="526" r="5" fill="#757575"/><text class="legend" x="311" y="530">Does not advance</text>','</svg>']
    path.write_text("\n".join(parts)+"\n",encoding="utf-8")
def _combine(root:Path,prepared:dict[str,Any],commands:dict[str,list[str]])->dict[str,Any]:
    predictions={};sources={}
    for identifier in prepared["record_ids"]:
        output=Path(commands[identifier][commands[identifier].index("--out_dir")+1]);source=_find_affinity(output,identifier);data=_load(source);prob=data.get("affinity_probability_binary");value=data.get("affinity_pred_value")
        if not isinstance(prob,(int,float)) or not math.isfinite(prob) or not 0<=prob<=1:raise ValueError(f"invalid affinity_probability_binary for {identifier}: {prob!r}")
        predictions[identifier]=(float(prob),float(value) if isinstance(value,(int,float)) and math.isfinite(value) else None);sources[identifier]=source
    reference=predictions["reference"][0];margin=float(prepared["comparison_margin"]);rows=[]
    for identifier in prepared["candidate_ids"]:
        probability,value=predictions[identifier];passed=probability>reference+margin;rows.append({"candidate_id":identifier,"affinity_probability_binary":probability,"affinity_pred_value_log10_ic50_uM":value,"reference_affinity_probability_binary":reference,"comparison_margin":margin,"decision":"advance_to_structural_confidence" if passed else "does_not_advance","detail":"strictly_greater_than_reference_plus_margin" if passed else "not_strictly_greater_than_reference_plus_margin","source_affinity_json":str(sources[identifier].relative_to(root))})
    dest=root/"boltz2_affinity/results";dest.mkdir(parents=True,exist_ok=True);table=dest/"boltz2_affinity.csv"
    with table.open("w",newline="",encoding="utf-8") as handle:w=csv.DictWriter(handle,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    qualified=[r for r in rows if r["decision"]=="advance_to_structural_confidence"];sequences={r.header.split()[0]:r.sequence for r in read_fasta(_candidate_fasta(root))};fasta=dest/"affinity_qualified_candidates.fasta";fasta.write_text("".join(f">{r['candidate_id']}\n{sequences[r['candidate_id']]}\n" for r in qualified),encoding="utf-8")
    plot=dest/"boltz2_affinity.svg";_write_plot(plot,rows,reference,margin)
    summary={"schema_version":1,"metric":"affinity_probability_binary","direction":"higher_is_better","strict_comparison":True,"reference_affinity_probability_binary":reference,"comparison_margin":margin,"candidate_count":len(rows),"qualified_count":len(qualified),"qualified_ids":[r["candidate_id"] for r in qualified],"plot":str(plot.relative_to(root)),"interpretation_warning":"Boltz-2 predictions do not establish experimental binding."};summary_path=dest/"summary.json";_write_json(summary_path,summary)
    return {**summary,"table":table,"fasta":fasta,"plot_path":plot,"summary_path":summary_path}
def run_boltz2_affinity(project_dir:Path,*,config_path:Path,prepare_only:bool=False)->int:
    root=project_dir.expanduser().resolve();prepared=prepare_boltz2_affinity(root);commands=build_boltz2_commands(root,config_path=config_path);start_stage(root,"boltz2_affinity",command=[x for identifier in prepared["record_ids"] for x in commands[identifier]]);log_path=root/"logs/boltz2_affinity.log";log_path.parent.mkdir(parents=True,exist_ok=True)
    try:
        if prepare_only:record_partial_stage_success(root,"boltz2_affinity",message="Boltz-2 YAML inputs prepared; inference not executed",outputs=[str(prepared["manifest_path"].relative_to(root))]);print("Preparation completed. Boltz-2 affinity remains resumable.");return 0
        with log_path.open("a",encoding="utf-8") as log:
            for identifier in prepared["record_ids"]:
                command=commands[identifier];output=Path(command[command.index("--out_dir")+1])
                if _completed(root,identifier,output):log.write(f"\nRESUME: reusing {identifier}\n");continue
                _mark(root,identifier,"running",command=command);log.write("\nCOMMAND: "+shlex.join(command)+"\n");log.flush();process=subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
                assert process.stdout is not None
                for line in process.stdout:print(line,end="");log.write(line);log.flush()
                process.stdout.close()
                if process.wait():raise RuntimeError(f"Boltz-2 failed for {identifier}")
                source=_find_affinity(output,identifier);_mark(root,identifier,"completed",command=command,affinity_json=str(source.relative_to(root)))
        result=_combine(root,prepared,commands);outputs=[str(log_path.relative_to(root)),str(result["table"].relative_to(root)),str(result["fasta"].relative_to(root)),str(result["plot_path"].relative_to(root)),str(result["summary_path"].relative_to(root))];finish_stage(root,"boltz2_affinity",succeeded=True,message=f"{result['qualified_count']} of {result['candidate_count']} candidates exceeded the reference affinity probability",outputs=outputs);print(f"Boltz-2 affinity completed: {result['qualified_count']} of {result['candidate_count']} candidates advance.")
        if result["qualified_count"]==0:
            from .structural_confidence import finalize_zero_survivor_workflow
            report=finalize_zero_survivor_workflow(root);print("No affinity-qualified candidates: structural confidence skipped.");print(f"Terminal report: {report['markdown']}")
        return 0
    except Exception as error:finish_stage(root,"boltz2_affinity",succeeded=False,message=str(error));raise
