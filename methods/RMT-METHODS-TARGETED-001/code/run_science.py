"""One-command server execution, immutable output directory and honest receipts."""
import argparse, hashlib, json, subprocess, sys
from datetime import datetime, timezone
from pathlib import Path

def main():
    a=argparse.ArgumentParser(); a.add_argument('--workspace',required=True,type=Path); a.add_argument('--output',required=True,type=Path); x=a.parse_args()
    root=Path(__file__).resolve().parents[1]; dest=x.output.resolve(); dest.mkdir(parents=True,exist_ok=False)
    def sha(p):
        with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
    before={'started_utc':datetime.now(timezone.utc).isoformat(),'python':sys.version,'files':{p.name:sha(p) for p in sorted((root/'code').glob('*.py'))},'protocol_sha256':sha(root/'PROTOCOL.json')}
    with (dest/'IMPLEMENTATION_FREEZE.json').open('x') as f:json.dump(before,f,indent=2)
    commands=[[sys.executable,'-m','unittest','test_analysis','-v'],[sys.executable,str(root/'code/analyze.py'),'--workspace',str(x.workspace),'--output',str(dest/'results')]]
    rc=0; receipts=[]
    for i,cmd in enumerate(commands):
        p=subprocess.run(cmd,cwd=root/'code',capture_output=True,text=True)
        with (dest/f'step_{i+1}.log').open('x') as f:f.write(p.stdout+p.stderr)
        receipts.append({'command':cmd,'exit_code':p.returncode,'log_sha256':sha(dest/f'step_{i+1}.log')})
        print(p.stdout+p.stderr,flush=True)
        if p.returncode:rc=p.returncode;break
    with (dest/'EXECUTION.json').open('x') as f:json.dump({'completed_utc':datetime.now(timezone.utc).isoformat(),'steps':receipts,'exit_code':rc},f,indent=2)
    raise SystemExit(rc)

if __name__=='__main__':main()
