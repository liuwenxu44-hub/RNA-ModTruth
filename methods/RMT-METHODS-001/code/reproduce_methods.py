"""Bounded offline replay from saved joins; no BAM/caller/GPU."""
import argparse, hashlib, json, os, subprocess, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
M=ROOT/'methods/RMT-METHODS-001'
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',required=True)
    p.add_argument('--check-only',action='store_true')
    a=p.parse_args()
    out=Path(a.output).resolve()
    if not out.is_relative_to(M/'results') or out==M/'results' or out.exists():
        raise SystemExit('Output must be a new absent child directory of methods/RMT-METHODS-001/results.')
    inputs=json.loads((M/'SERVER_INPUT_MANIFEST.json').read_text())
    for r in inputs['files']:
        q=ROOT/r['path']
        if not q.is_file() or sha(q)!=r['sha256']:raise SystemExit('Missing or changed input: '+r['path'])
    commands=[['code/test_methods.py'],['code/robustness.py','--output',str(out/'robustness')],
              ['code/check_regression.py','--output',str(out/'facts')]]
    if a.check_only:
        print(json.dumps(dict(status='INPUT_AND_REPLAY_PLAN_VERIFIED_NOT_EXECUTED',input_files=len(inputs['files']),commands=commands),indent=2));return
    out.mkdir(parents=True,exist_ok=False)
    start=time.monotonic();env=os.environ.copy()
    env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    for i,cmd in enumerate(commands,1):
        with (out/f'step_{i}.log').open('x') as f:
            rc=subprocess.run([sys.executable,'-B',*cmd],cwd=M,env=env,stdout=f,stderr=subprocess.STDOUT).returncode
        if rc:raise SystemExit(f'Step {i} failed; inspect the new log. Historical results remain unchanged.')
    compared=[]
    for name in ['METRICS.json','METRICS.tsv','SAMPLING_FRAME.tsv','SELECTED_READS_LOCAL_ONLY.tsv','NATIVE_COUNTS.json','LEGACY_RECONCILIATION.json']:
        assert sha(out/'robustness'/name)==sha(M/'results/server_run_001/robustness'/name),name
        compared.append(name)
    for name in ['FACTS.json','METHOD_CHECKS.json','READ_MAP_REGRESSION.json']:
        assert json.loads((out/'facts'/name).read_text())==json.loads((M/'results/server_run_001/facts'/name).read_text()),name
        compared.append(name)
    report=dict(status='PASS',wall_seconds=time.monotonic()-start,threads=1,compared_files=compared,
                no_network=True,no_caller=True,no_BAM_reprocessing=True,new_operating_system_validation=False)
    with (out/'REPLAY.json').open('x') as f:json.dump(report,f,indent=2)
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
