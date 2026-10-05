"""Transparent downstream consumer, independent of every validation decision.

It deliberately executes the supplied join/population/point-score contract so an
injected error can produce biased values, missing joins or no result. Those are
diagnostic counterfactuals, never accepted biological estimates.
"""
import argparse
import collections
import hashlib
import json
import math
from bundle import load, write_json


def calculate(path):
    c, rows, truth, refs, mapping=load(path)
    keys=c['truth_join_keys']
    labels=collections.defaultdict(list)
    for r in truth:
        normalized=dict(r,position=str(int(r['position'])-c['truth_origin']))
        labels[tuple(normalized[k] for k in keys)].append(r)
    counts=collections.Counter()
    counts['input_records']=len(rows)
    counts['input_unique_evidence_ids']=len({r['evidence_id'] for r in rows})
    counts['declared_population_opportunities']=c['expected_opportunities']
    confusion=collections.Counter(dict(TP=0,TN=0,FP=0,FN=0))
    scored=[]
    normalized_records=[]
    for r in rows:
        coordinate=str(int(r['position'])-c['consumer_origin'])
        normalized=dict(r,position=coordinate)
        match=labels.get(tuple(normalized[k] for k in keys),[])
        counts['truth_join_candidates']+=len(match)
        if len(match)==0:
            counts['no_truth_join']+=1
            continue
        if len(match)!=1:
            counts['ambiguous_truth_join']+=1
            continue
        counts['unique_truth_join_records']+=1
        t=match[0]
        normalized_records.append((r['evidence_id'],r['sample_id'],r['run_id'],r['read_id'],
            c['contig_aliases'].get(r['contig'],r['contig']),coordinate,r['strand'],
            r['modification_code'],r['score_kind'],r['score_status'],r['p_mid'],t['truth_label']))
        if r['alignment_status']!='MATCHED_READ_SITE':
            counts['alignment_ineligible']+=1
            continue
        counts['eligible_aligned_A_sites']+=1
        if not r['p_mid']:
            counts['eligible_missing_probability']+=1
            continue
        p=float(r['p_mid'])
        y=int(t['truth_label'])
        pred=int(p>=.5)
        counts['evaluable_points']+=1
        counts['evaluable_control' if y==0 else 'evaluable_modified']+=1
        confusion['TP' if y and pred else 'FN' if y else 'FP' if pred else 'TN']+=1
        lo=float(r['p_low']) if r['p_low'] else None
        hi=float(r['p_high']) if r['p_high'] else None
        # The interval is applicable only if the supplied point is its declared
        # midpoint and the semantics remain read-site MM/ML, without renormalization.
        interval_ok=(lo is not None and hi is not None and r['ml_integer']!=''
            and r['score_kind']=='read_site_ML_interval' and r['modification_code'] in ('a','28871')
            and p==(int(r['ml_integer'])+.5)/256 and lo==int(r['ml_integer'])/256 and hi==(int(r['ml_integer'])+1)/256)
        if not interval_ok:
            counts['points_without_valid_native_interval']+=1
        scored.append((p,y,lo,hi,interval_ok))
    def loss(p,y):
        q=p if y else 1-p
        return -math.log(q) if q>0 else math.inf
    n=len(scored)
    brier=sum((p-y)**2 for p,y,_,_,_ in scored)/n if n else None
    logloss=sum(loss(p,y) for p,y,_,_,_ in scored)/n if n else None
    fpr=confusion['FP']/(confusion['FP']+confusion['TN']) if confusion['FP']+confusion['TN'] else None
    tpr=confusion['TP']/(confusion['TP']+confusion['FN']) if confusion['TP']+confusion['FN'] else None
    intervals=None
    if n and all(x[4] for x in scored):
        bmin=sum(min((lo-y)**2,(hi-y)**2) for _,y,lo,hi,_ in scored)/n
        bmax=sum(max((lo-y)**2,(hi-y)**2) for _,y,lo,hi,_ in scored)/n
        llmin=sum(min(loss(lo,y),loss(hi,y)) for _,y,lo,hi,_ in scored)/n
        llmax=sum(max(loss(lo,y),loss(hi,y)) for _,y,lo,hi,_ in scored)/n
        # Upper endpoint is excluded: touching 0.5 from below does not cross it.
        ambiguous=sum(lo<.5<hi for _,_,lo,hi,_ in scored)
        intervals=dict(Brier_infimum=bmin,Brier_supremum=bmax,
                       log_loss_infimum=llmin,log_loss_supremum=llmax if math.isfinite(llmax) else 'INFINITE_SUPREMUM',
                       threshold_interval_ambiguous=ambiguous,
                       meaning='bounds within original ML intervals; not sampling confidence intervals')
    for k in ('truth_join_candidates','no_truth_join','ambiguous_truth_join','unique_truth_join_records',
              'eligible_aligned_A_sites','eligible_missing_probability','alignment_ineligible','evaluable_points',
              'evaluable_control','evaluable_modified','points_without_valid_native_interval'):
        counts.setdefault(k,0)
    return dict(counts=dict(counts),confusion=dict(confusion),fixed_threshold=.5,
        balanced_accuracy=(tpr+1-fpr)/2 if fpr is not None and tpr is not None else None,FPR=fpr,
        Brier=brier,log_loss=logloss if logloss is None or math.isfinite(logloss) else 'INFINITE',
        valid_probability_semantics=not counts['points_without_valid_native_interval'],
        probability_metrics_status='SCOPED_NATIVE_INTERVAL' if n and not counts['points_without_valid_native_interval'] else 'NO_EVALUABLE_RESULT' if not n else 'INVALID_INPUT_DIAGNOSTIC_ONLY',
        ML_interval_sensitivity=intervals,
        explicit_coverage_of_claimed_opportunities=n/c['expected_opportunities'] if c['expected_opportunities'] else None,
        filter_coverage_records=len(rows)/c['expected_opportunities'] if c['expected_opportunities'] else None,
        truth_join_yield=counts['unique_truth_join_records']/len(rows) if rows else None,
        normalized_source_join_sha256=hashlib.sha256(json.dumps(sorted(normalized_records),separators=(',',':')).encode()).hexdigest(),
        scope='within this provider-supplied selected construct bundle; not native occupancy or between-model comparison')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--bundle',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    write_json(a.output,calculate(a.bundle))
