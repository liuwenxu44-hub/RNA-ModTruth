"""New, symmetric read-map regression tests; historical versions stay read-only."""
import argparse, csv, importlib.util, json, math, sys, time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]; HERE=Path(__file__).resolve().parents[1]; E=ROOT/'execution/RMT-EXEC-001'
sys.path.insert(0,str(E/'code'))
from bundle import load,write_json,write_table,digest
from metrics import calculate
def module(name,path):
    sp=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(sp);sp.loader.exec_module(m);return m
newr=module('newr',HERE/'code/checkers_v2/rmt.py'); newb=module('newb',HERE/'code/checkers_v2/simple_baseline.py')
oldr=module('oldr',E/'code/rmt.py'); oldb=module('oldb',E/'code/simple_baseline.py')

ROLES=[('DEVELOPMENT',E/'canonical_bundles/development_v2',E/'controlled_cases/development_v2',E/'results/development_test_final'),
       ('EVALUATION',E/'results/evaluation_run_v1/canonical_bundles/evaluation',E/'results/evaluation_run_v1/controlled_cases/evaluation',E/'results/evaluation_run_v1/results/evaluation')]
def equal(a,b):
    if isinstance(a,dict): return a.keys()==b.keys() and all(equal(a[k],b[k]) for k in a)
    if isinstance(a,list):return len(a)==len(b) and all(equal(x,y) for x,y in zip(a,b))
    if isinstance(a,float) and isinstance(b,(int,float)):return math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)
    return a==b
def subclass_metrics(bundle):
    c,rows,truth,_,_=load(bundle); keys=c['truth_join_keys']; index={}
    for t in truth:
        x=dict(t,position=str(int(t['position'])-c['truth_origin'])); index.setdefault(tuple(x[k] for k in keys),[]).append(t)
    n=[0,0]; s=[0.,0.]; counts={str(y):dict(truth_overlap=0,eligible=0,native_explicit=0,scored=0) for y in (0,1)}
    for r in rows:
        x=dict(r,position=str(int(r['position'])-c['consumer_origin'])); ts=index.get(tuple(x[k] for k in keys),[])
        if len(ts)!=1:continue
        y=int(ts[0]['truth_label']); z=counts[str(y)]; z['truth_overlap']+=1
        if r['alignment_status']!='MATCHED_READ_SITE':continue
        z['eligible']+=1; z['native_explicit']+=r['ml_integer']!=''
        if not r['p_mid']:continue
        n[y]+=1; s[y]+=(float(r['p_mid'])-y)**2;z['scored']+=1
    b=[s[y]/n[y] if n[y] else None for y in (0,1)]
    return dict(class_denominators=counts,Brier_control=b[0],Brier_modified=b[1],Brier_equal_class=sum(b)/2 if all(x is not None for x in b) else None)

def run(out):
    start=time.monotonic();out.mkdir(parents=True,exist_ok=False);checks=[];facts=[];reg=[]
    for role,parent,cases,results in ROLES:
        cases_index=json.loads((cases/'CASE_INDEX.json').read_text())
        with (results/'method_results.tsv').open() as f: methodrows=list(csv.DictReader(f,delimiter='\t'))
        for case in cases_index:
            cid=case['case_id']; actual=calculate(cases/cid); frozen=json.loads((results/'metrics'/f'{cid}.json').read_text())
            assert equal(actual,frozen),(role,cid,'metric drift')
            fact=dict(role=role,case_id=cid,evidence_class={'CANONICAL':'OBSERVED_NATIVE_INPUT','CONTROLLED_ERROR':'CONTROLLED_ERROR_INJECTION','LEGAL_TRANSFORM':'LEGAL_TRANSFORMATION'}[case['kind']],source_bundle=str((cases/cid).relative_to(ROOT)),source_metric=str((results/'metrics'/f'{cid}.json').relative_to(ROOT)),source_metric_sha256=digest(results/'metrics'/f'{cid}.json'),recalculated_metrics=actual,**subclass_metrics(cases/cid))
            facts.append(fact)
            for method in ('B0','B1','B2','RMT'):
                raw=json.loads((results/'checks'/f'{cid}_{method}.json').read_text()); tab=next(r for r in methodrows if r['case_id']==cid and r['method']==method)
                assert raw['status']==tab['actual_status'],(role,cid,method)
                check=dict(role=role,case_id=cid,method=method,status=raw['status'],reported_detects_error=tab['detects_error'],reported_legal_false_rejection=tab['legal_false_rejection'],recovery=raw.get('recovery',{}),source_check=str((results/'checks'/f'{cid}_{method}.json').relative_to(ROOT)),source_check_sha256=digest(results/'checks'/f'{cid}_{method}.json'))
                if method in ('B2','RMT'):
                    calc=(newb.check if method=='B2' else newr.validate)(cases/cid,parent)
                    assert calc['status']==raw['status'] and calc['issues']==raw['issues'],(role,cid,method,'new regression changes original outcomes')
                    check['new_v2_status']=calc['status'];check['new_v2_issue_identity']='UNCHANGED'
                    if raw.get('recovery',{}).get('written'):
                        rec=results/'recovered'/f'{cid}_{method}'; rm=calculate(rec)
                        oldm=json.loads((results/'metrics'/f'{cid}_{method}_recovered.json').read_text())
                        assert equal(rm,oldm) and equal(rm,calculate(parent))
                        facts.append(dict(role=role,case_id=cid+'_'+method+'_recovered',evidence_class='SAFE_RECOVERY',source_bundle=str(rec.relative_to(ROOT)),source_metric=str((results/'metrics'/f'{cid}_{method}_recovered.json').relative_to(ROOT)),source_metric_sha256=digest(results/'metrics'/f'{cid}_{method}_recovered.json'),recalculated_metrics=rm,**subclass_metrics(rec)))
                checks.append(check)
        c,_,_,_,mapping=load(parent)
        for test in ('canonical','read_map_reverse','duplicate_identical_key','conflicting_key','missing_row','sample_identity_change','payload_change'):
            m=[dict(x) for x in mapping]
            if test=='read_map_reverse':m.reverse()
            elif test=='duplicate_identical_key':m.append(dict(m[0]))
            elif test=='conflicting_key':
                x=dict(m[0]);x['read_id']='CONTROLLED_CONFLICT';m.append(x)
            elif test=='missing_row':m=m[1:]
            elif test=='sample_identity_change':m[0]['sample_id']='CONTROLLED_WRONG_SAMPLE'
            elif test=='payload_change':m[0]['source_bam_sha256']='0'*64
            folder=out/'fixtures'/role/test;folder.mkdir(parents=True,exist_ok=False)
            for name in c['prediction_parts']+['truth.tsv','reference.fa','contract.json']:(folder/name).symlink_to(parent/name)
            write_table(folder/'read_map.tsv',m)
            for method,newfun,oldfun in [('B2',newb.check,oldb.check),('RMT',newr.validate,oldr.validate)]:
                new=newfun(folder,parent); expected='ACCEPT' if test in ('canonical','read_map_reverse') else 'UNDETERMINED'
                old=oldfun(folder,parent)['status'] if test in ('canonical','read_map_reverse') else 'NOT_RERUN'
                assert new['status']==expected,(role,test,method,new['status'])
                reg.append(dict(role=role,test=test,method=method,old_frozen_status=old,new_v2_status=new['status'],expected=expected,status='PASS',evidence_class='LEGAL_TRANSFORMATION' if test in ('canonical','read_map_reverse') else 'CONTROLLED_ERROR_INJECTION',analysis_role='NEW_SOFTWARE_REGRESSION_NOT_ORIGINAL_EVALUATION'))
        print(role,'historical facts and symmetric regression verified',flush=True)
    write_json(out/'FACTS.json',facts);write_json(out/'METHOD_CHECKS.json',checks);write_json(out/'READ_MAP_REGRESSION.json',reg)
    write_json(out/'AUDIT.json',dict(status='PASS',historical_metric_checks=len(facts),method_source_checks=len(checks),regression_checks=len(reg),wall_seconds=time.monotonic()-start,old_frozen_code_changed=False,B2_RMT_support_equal=True))
    print('PASS',len(facts),len(checks),len(reg),flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',required=True);a=p.parse_args();run(Path(a.output))
