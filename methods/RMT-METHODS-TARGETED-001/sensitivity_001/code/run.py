"""Single bounded test and analysis invocation with immutable receipts."""
import argparse,hashlib,json,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path

p=argparse.ArgumentParser();p.add_argument('--workspace',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
r=Path(__file__).resolve().parents[1];out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
def sha(f):
    with f.open('rb') as x:return hashlib.file_digest(x,'sha256').hexdigest()
freeze={'started_utc':datetime.now(timezone.utc).isoformat(),'python':sys.version,'protocol_sha256':sha(r/'PROTOCOL.json'),'files':{name:sha(r/'code'/name) for name in ['analyze_intervals.py','test_intervals.py','run.py']}}
with (out/'IMPLEMENTATION_FREEZE.json').open('x') as x:json.dump(freeze,x,indent=2)
cmds=[[sys.executable,'-m','unittest','test_intervals','-v'],[sys.executable,str(r/'code/analyze_intervals.py'),'--workspace',str(a.workspace),'--output',str(out/'results')]];steps=[];rc=0
for i,cmd in enumerate(cmds,1):
    z=subprocess.run(cmd,cwd=r/'code',capture_output=True,text=True)
    log=out/f'step_{i}.log'
    with log.open('x') as x:x.write(z.stdout+z.stderr)
    steps.append({'command':cmd,'exit_code':z.returncode,'log_sha256':sha(log)});print(z.stdout+z.stderr,flush=True)
    if z.returncode:rc=z.returncode;break
with (out/'EXECUTION.json').open('x') as x:json.dump({'completed_utc':datetime.now(timezone.utc).isoformat(),'steps':steps,'exit_code':rc},x,indent=2)
raise SystemExit(rc)
