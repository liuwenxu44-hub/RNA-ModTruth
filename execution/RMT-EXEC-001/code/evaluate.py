"""Batch evaluation; no stop on the first failing method/case; no checker mutation."""
import argparse
import collections
import json
import math
import sys
from pathlib import Path
from bundle import read_json,write_json,write_table,digest
from measure import measured
from metrics import calculate

METHODS={'B0':'structure_baseline.py','B1':'tools_baseline.py','B2':'simple_baseline.py','RMT':'rmt.py'}
B1_UNASSESSED={'site_score_as_read_probability','filtered_output_claimed_all_records','missing_probability_imputed_zero'}


def equivalent(a,b):
    if isinstance(a,dict) and isinstance(b,dict):
        return set(a)==set(b) and all(equivalent(a[k],b[k]) for k in a)
    if isinstance(a,list) and isinstance(b,list):
        return len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if type(a) is float and type(b) is float:
        return math.isclose(a,b,abs_tol=1e-12,rel_tol=1e-12)
    return a==b


def evaluate(cases,evidence,output):
    cases,evidence,output=Path(cases).resolve(),Path(evidence).resolve(),Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    (output/'checks').mkdir()
    (output/'logs').mkdir()
    (output/'recovered').mkdir()
    (output/'metrics').mkdir()
    index=read_json(cases/'CASE_INDEX.json')
    canonical=calculate(cases/'00_canonical')
    records=[]
    runtimes=[]
    downstream=[]
    for case in index:
        cid=case['case_id']
        path=cases/cid
        for file,sha in case['generated_files_sha256'].items():
            if digest(path/file)!=sha:
                raise ValueError('CASE_HASH_DRIFT:'+cid+':'+file)
        raw=calculate(path)
        write_json(output/'metrics'/f'{cid}.json',raw)
        downstream.append(dict(case_id=cid,kind=case['kind'],method='UNVALIDATED_CONSUMER',metrics=raw,
            legal_semantic_invariance=equivalent(raw,canonical) if case['kind']=='LEGAL_TRANSFORM' else None))
        for method,script in METHODS.items():
            target=output/'checks'/f'{cid}_{method}.json'
            command=[sys.executable,'-B',Path(__file__).with_name(script),'--bundle',path,'--evidence',evidence,'--output',target]
            recovered=output/'recovered'/f'{cid}_{method}'
            if method in ('B2','RMT'):
                command+=['--recover',recovered]
            timing=measured(command,output/'logs'/f'{cid}_{method}.log')
            timing.update(case_id=cid,method=method)
            runtimes.append(timing)
            if timing['exit_code'] or not target.exists():
                result=dict(status='EXECUTION_ERROR',issues=[],recovery={'status':'NOT_EXECUTED','written':False})
            else:
                result=read_json(target)
            expected=set(case['affected_evidence_ids'])
            located=set(eid for item in result['issues'] for eid in item['evidence_ids'])
            supported=method!='B0' and not (case['kind']=='CONTROLLED_ERROR' and method=='B1' and case['family'] in B1_UNASSESSED)
            recovery=result.get('recovery',dict(status='UNSUPPORTED',written=False))
            recovery_equivalence=None
            if recovery.get('written'):
                fixed=calculate(recovered)
                write_json(output/'metrics'/f'{cid}_{method}_recovered.json',fixed)
                downstream.append(dict(case_id=cid,kind='SAFE_RECOVERY',method=method,metrics=fixed))
                recovery_equivalence=equivalent(fixed,canonical)
            row=dict(case_id=cid,family=case['family'],kind=case['kind'],role=case['role'],method=method,
                actual_status=result['status'],endpoint_supported=supported,
                endpoint_result=('SUPPORTED_SCOPE_ACCEPT' if method=='B1' and result['status']=='ACCEPT' else result['status']) if supported else 'UNASSESSED',
                detects_error=(result['status']=='REJECT') if case['kind']=='CONTROLLED_ERROR' and supported else None,
                legal_false_rejection=(result['status']!='ACCEPT') if case['kind']!='CONTROLLED_ERROR' else None,
                affected_records=len(expected),localized_records=len(located),
                localization_recall=len(located&expected)/len(expected) if expected and case['kind']=='CONTROLLED_ERROR' else None,
                localization_precision=len(located&expected)/len(located) if located and expected and case['kind']=='CONTROLLED_ERROR' else None,
                all_explanations_have_locators=all(i.get('locator') and i.get('evidence_locator') for i in result['issues']) if result['issues'] else None,
                recovery_status=recovery['status'],recovery_equivalent=recovery_equivalence,
                wall_seconds=timing['wall_seconds'],peak_rss_kib=timing['peak_rss_kib'])
            records.append(row)
        print(json.dumps({'case':cid,'role':case['role'],'completed_methods':4}),flush=True)
    write_table(output/'method_results.tsv',records)
    write_json(output/'runtime_detail.json',runtimes)
    write_json(output/'downstream_results.json',downstream)
    summary={}
    for method in METHODS:
        r=[x for x in records if x['method']==method]
        err=[x for x in r if x['kind']=='CONTROLLED_ERROR']
        legal=[x for x in r if x['kind']=='LEGAL_TRANSFORM']
        summary[method]=dict(error_cases=len(err),supported_error_cases=sum(x['endpoint_supported'] for x in err),
            detected=sum(x['detects_error'] is True for x in err),unsupported_errors=sum(not x['endpoint_supported'] for x in err),
            legal_cases=len(legal),legal_false_rejections=sum(x['legal_false_rejection'] for x in legal),
            undetermined=sum(x['actual_status']=='UNDETERMINED' for x in r),execution_errors=sum(x['actual_status']=='EXECUTION_ERROR' for x in r),
            recoveries=sum(x['recovery_equivalent'] is not None for x in r),non_equivalent_recoveries=sum(x['recovery_equivalent'] is False for x in r),
            total_wall_seconds=sum(x['wall_seconds'] for x in r),peak_rss_kib=max(x['peak_rss_kib'] for x in r))
    write_json(output/'SUMMARY.json',summary)
    print(json.dumps(summary),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('cases','evidence','output'):
        p.add_argument('--'+k,required=True)
    a=p.parse_args()
    evaluate(a.cases,a.evidence,a.output)
