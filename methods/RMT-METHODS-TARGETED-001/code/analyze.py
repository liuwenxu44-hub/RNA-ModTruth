"""Bounded analysis of frozen saved joins; stdlib, one streaming read per file."""
import argparse, collections, csv, gzip, hashlib, json, math, platform, resource, time
from datetime import datetime, timezone
from pathlib import Path

TOL=1e-10
def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def write_json(p,x):
    with p.open('x') as f: json.dump(x,f,indent=2,allow_nan=False)
def write_tsv(p,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('x',newline='') as f:
        w=csv.DictWriter(f,keys,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def ratio(a,b): return a/b if b else None
def average(vals): return sum(vals)/len(vals) if vals and all(v is not None for v in vals) else None

def retention(m,policy,t):
    lo,hi,p=m/256,(m+1)/256,(m+.5)/256
    if policy=='explicit_baseline': return True,True,True
    if policy=='one_sided_high_p_reference': return p>=t,lo>=t,hi>t
    if policy!='symmetric_confidence': raise ValueError(policy)
    if t==.5: return True,True,True
    return max(p,1-p)>=t,(lo>=t or hi<=1-t),(hi>t or lo<=1-t)

class Group:
    def __init__(self): self.opp=0; self.eligible=0; self.hist=collections.Counter(); self.reasons=collections.Counter()
    def metric(self,y,policy,t):
        n=b=err=definite=possible=0
        for m,k in self.hist.items():
            keep,d,z=retention(m,policy,t); definite+=d*k; possible+=z*k
            if keep:
                p=(m+.5)/256; n+=k; b+=k*(p-y)**2; err+=k*(int(p>=.5)!=y)
        o=sum(self.hist.values())
        return dict(opportunities=self.opp,N=self.eligible,O=o,S=n,coverage_S_over_O=ratio(n,o),coverage_S_over_N=ratio(n,self.eligible),brier=ratio(b,n),risk=ratio(err,n),loss_sum=b,error_count=err,definite_retained=definite,possible_retained=possible,ambiguous_retention=possible-definite,state='ESTIMABLE' if n else 'NO_RETAINED_RECORDS')

def pooled(classrows):
    out={}
    for key in ['opportunities','N','O','S','loss_sum','error_count','definite_retained','possible_retained','ambiguous_retention']:
        out[key]=sum(r[key] for r in classrows)
    b=[r['brier'] for r in classrows]; e=[r['risk'] for r in classrows]
    out.update(Brier_observed=ratio(out['loss_sum'],out['S']),Brier_equal_class=average(b),risk_observed=ratio(out['error_count'],out['S']),risk_equal_class=average(e),balanced_accuracy=None if average(e) is None else 1-average(e),coverage_S_over_O=ratio(out['S'],out['O']),coverage_S_over_N=ratio(out['S'],out['N']),state='ESTIMABLE' if all(r['S'] for r in classrows) else 'CLASS_MISSING')
    for y,r in enumerate(classrows):
        for k in ['N','O','S','brier','risk','coverage_S_over_O','coverage_S_over_N','ambiguous_retention']: out[f'{k}_class_{y}']=r[k]
    return out

def decompose(base,new,metadata):
    """Inputs keyed strata with (scored n, mean loss); no missing-loss imputation."""
    allkeys=sorted(base); support=[k for k in allkeys if base[k][0] and new[k][0]]
    n0=sum(v[0] for v in base.values()); n1=sum(v[0] for v in new.values())
    full0=sum(n*l for n,l in base.values() if n)/n0 if n0 else None
    full1=sum(n*l for n,l in new.values() if n)/n1 if n1 else None
    missing=[k for k in allkeys if k not in support]
    out=[]
    scopes=[('FULL_SUPPORT',allkeys)] if not missing else [('FULL_SUPPORT',None),('COMMON_SUPPORT_ONLY',support)]
    for scope,keys in scopes:
        r={**metadata,'scope':scope,'excluded_strata':'|'.join(missing),'excluded_baseline_weight':ratio(sum(base[k][0] for k in missing),n0),'excluded_selected_weight':ratio(sum(new[k][0] for k in missing),n1),'full_observed_delta':None if full0 is None or full1 is None else full1-full0}
        if not keys:
            r.update(state='UNIDENTIFIABLE',baseline_loss=None,selected_loss=None,delta=None,within_loss_component=None,weight_component=None,addback_residual=None)
        else:
            a=sum(base[k][0] for k in keys); b=sum(new[k][0] for k in keys)
            terms=[]
            for k in keys:
                w0,w1=base[k][0]/a,new[k][0]/b; l0,l1=base[k][1],new[k][1]
                terms.append(((w1+w0)/2*(l1-l0),(l1+l0)/2*(w1-w0)))
            loss0=sum(base[k][0]*base[k][1] for k in keys)/a; loss1=sum(new[k][0]*new[k][1] for k in keys)/b
            within=sum(x[0] for x in terms); weight=sum(x[1] for x in terms); residual=loss1-loss0-within-weight
            assert abs(residual)<TOL
            r.update(state='IDENTIFIED',baseline_loss=loss0,selected_loss=loss1,delta=loss1-loss0,within_loss_component=within,weight_component=weight,addback_residual=residual)
        out.append(r)
    return out

def hiding(g,y):
    n=sum(g.hist.values()); hidden=0; low=high=truth=0.
    for m,k in g.hist.items():
        lo,hi,p=m/256,(m+1)/256,(m+.5)/256
        if hi<=.25: lo,hi=0.,.25; hidden+=k
        low+=k*(lo**2 if y==0 else (1-hi)**2)
        high+=k*(hi**2 if y==0 else (1-lo)**2)
        truth+=k*(p-y)**2
    assert not n or low-TOL<=truth<=high+TOL
    return dict(n=n,hidden=hidden,lower_sum=low,upper_sum=high,known_loss_sum=truth)

def run(workspace,out):
    root=Path(__file__).resolve().parents[1]; freeze=json.loads((root/'INPUTS.json').read_text()); protocol=json.loads((root/'PROTOCOL.json').read_text())
    assert sha(root/'PROTOCOL.json')==freeze['protocol_sha256']
    out.mkdir(exist_ok=False,parents=True); started=time.monotonic(); cpu=time.process_time()
    groups=collections.defaultdict(Group); sites=collections.defaultdict(lambda: collections.Counter()); counts=[]
    for inp in freeze['record_inputs']:
        p=workspace/inp['path']; assert sha(p)==inp['sha256']
        release='2025_DRACH' if inp['role'].startswith('development') else '2026_all5mer'
        n=0
        with gzip.open(p,'rt',newline='') as f:
            for r in csv.DictReader(f,delimiter='\t'):
                n+=1; y=int(r['truth_label']); assert y in (0,1)
                g=groups[(release,r['contig'],y)]; g.opp+=1; g.reasons[(r['alignment_status'],r['score_status'])]+=1
                eligible=r['alignment_status']=='MATCHED_READ_SITE'; g.eligible+=eligible
                if not r['p_mid']: assert r['score_status']!='explicit'; continue
                m=int(r['ml_integer']); assert 0<=m<=255 and eligible and r['score_status']=='explicit'
                assert r['score_kind']=='read_site_ML_interval' and r['modification_code'] in ('a','28871')
                assert float(r['p_mid'])==(m+.5)/256 and float(r['p_low'])==m/256 and float(r['p_high'])==(m+1)/256
                g.hist[m]+=1
                sk=(release,r['sample_id'],r['run_id'],r['contig'],r['ref_pos0'],r['strand'])
                sites[sk][(y,m)]+=1
        counts.append(dict(release=release,role=inp['role'],rows=n,sha256=inp['sha256'])); print(inp['role'],n,flush=True)
    classrows=[]; curve=[]; fixed=[]; decompositions=[]; identity=[]; boundrows=[]; reasonrows=[]; reconciliation=[]
    policies=[('explicit_baseline',.5)]+[('symmetric_confidence',t) for t in protocol['A']['thresholds']]+[('one_sided_high_p_reference',.9)]
    for release in sorted({k[0] for k in groups}):
        constructs=sorted({k[1] for k in groups if k[0]==release}); metrics={}; allgroups=[]
        for y in (0,1):
            a=Group()
            for c in constructs:
                g=groups[(release,c,y)]; a.opp+=g.opp; a.eligible+=g.eligible; a.hist.update(g.hist); a.reasons.update(g.reasons)
                for (alignment,score),n in sorted(g.reasons.items()): reasonrows.append(dict(release=release,construct=c,truth_class=y,alignment_status=alignment,score_status=score,count=n,eligible_for_native_bounds=alignment=='MATCHED_READ_SITE' and score!='explicit'))
            allgroups.append(a)
        for policy,t in policies:
            meta=dict(release=release,policy=policy,threshold=t,analysis_class=protocol['analysis_class'])
            for c in constructs+['ALL']:
                rr=[]
                for y in (0,1):
                    g=allgroups[y] if c=='ALL' else groups[(release,c,y)]; r=g.metric(y,policy,t)
                    classrows.append({**meta,'construct':c,'truth_class':y,**r}); rr.append(r)
                comb=pooled(rr); metrics[(policy,t,c)]=(rr,comb)
                curve.append({**meta,'construct':c,**comb})
            missing=[f'{c}:{y}' for c in constructs for y in (0,1) if not metrics[(policy,t,c)][0][y]['S']]
            fr={**meta,'K':len(constructs),'weight_definition':'1/K construct and 1/2 class within construct','missing_construct_class':'|'.join(missing),'state':'NOT_ESTIMABLE_FULL_SUPPORT' if missing else 'ESTIMABLE'}
            for key in ['Brier_equal_class','risk_equal_class']:
                fr['fixed_construct_'+key]=average([metrics[(policy,t,c)][1][key] for c in constructs])
            fixed.append(fr)
            if policy!='explicit_baseline':
                b0=metrics[('explicit_baseline',.5,'ALL')][0]; b1=metrics[(policy,t,'ALL')][0]
                decompositions+=decompose({str(y):(r['S'],r['brier']) for y,r in enumerate(b0)},{str(y):(r['S'],r['brier']) for y,r in enumerate(b1)},{**meta,'stratum_scheme':'truth_class'})
                b0={c:(metrics[('explicit_baseline',.5,c)][1]['S'],metrics[('explicit_baseline',.5,c)][1]['Brier_observed']) for c in constructs}
                b1={c:(metrics[(policy,t,c)][1]['S'],metrics[(policy,t,c)][1]['Brier_observed']) for c in constructs}
                decompositions+=decompose(b0,b1,{**meta,'stratum_scheme':'construct'})
        # Single controlled hiding demonstration; classes and fixed O weights.
        hh=[hiding(allgroups[y],y) for y in (0,1)]
        for label,h in [(str(y),v) for y,v in enumerate(hh)]+[('ALL',{k:sum(h[k] for h in hh) for k in hh[0]})]:
            boundrows.append(dict(release=release,scope='CONTROLLED_SCORE_REMOVAL',truth_class=label,denominator='original explicit O',N=h['n'],hidden=h['hidden'],loss_infimum=ratio(h['lower_sum'],h['n']),known_midpoint_loss=ratio(h['known_loss_sum'],h['n']),loss_supremum=ratio(h['upper_sum'],h['n']),contained=True))
        boundrows.append(dict(release=release,scope='CONTROLLED_SCORE_REMOVAL',truth_class='EQUAL_CLASS',denominator='1/2 per original explicit class',N=sum(h['n'] for h in hh),hidden=sum(h['hidden'] for h in hh),loss_infimum=average([ratio(h['lower_sum'],h['n']) for h in hh]),known_midpoint_loss=average([ratio(h['known_loss_sum'],h['n']) for h in hh]),loss_supremum=average([ratio(h['upper_sum'],h['n']) for h in hh]),contained=all(h['n'] for h in hh)>0))
        # Algebra identity tested at its exact grouping, never include label in group key.
        for c in constructs+['ALL']:
            selected=[h for k,h in sites.items() if k[0]==release and (c=='ALL' or k[3]==c)]
            n=mixed=0; before=after=variance=0.
            for h in selected:
                ys={y for y,m in h}; mixed+=len(ys)!=1
                count=sum(h.values()); mean=sum(k*(m+.5)/256 for (y,m),k in h.items())/count
                n+=count; before+=sum(k*((m+.5)/256-y)**2 for (y,m),k in h.items())
                after+=sum(k*(mean-y)**2 for (y,m),k in h.items())
                variance+=sum(k*((m+.5)/256-mean)**2 for (y,m),k in h.items())
            residual=ratio(before-after-variance,n) if not mixed else None
            if residual is not None: assert abs(residual)<TOL
            identity.append(dict(release=release,construct=c,explicit_records=n,site_groups=len(selected),mixed_label_groups=mixed,read_Brier=ratio(before,n),repeated_site_mean_Brier=ratio(after,n),population_within_site_variance=ratio(variance,n),identity_residual=residual,state='VERIFIED_UNIFORM_LABEL_IDENTITY' if not mixed else 'NOT_APPLICABLE_MIXED_LABEL'))
        parent=workspace/'methods/RMT-METHODS-001/results/server_run_001/robustness/METRICS.tsv'
        with parent.open() as f:
            old=list(csv.DictReader(f,delimiter='\t'))
        for condition,policy,t in [('canonical','explicit_baseline',.5),('high_p_filter','one_sided_high_p_reference',.9)]:
            for c in constructs+['ALL']:
                prior=next(r for r in old if r['release']==release and r['subset']=='full' and r['condition']==condition and r['construct']==c)
                curr=metrics[(policy,t,c)][1]
                for key in ['Brier_observed','Brier_equal_class']:
                    x=curr[key]; y=float(prior[key]) if prior[key] else None
                    ok=(x is None and y is None) or (x is not None and y is not None and abs(x-y)<TOL)
                    reconciliation.append(dict(release=release,construct=c,condition=condition,metric=key,old=y,new=x,match=ok)); assert ok,(release,c,condition,key,x,y)
        siteold=next(r for r in old if r['release']==release and r['subset']=='full' and r['condition']=='site_mean_as_read' and r['construct']=='ALL')
        current=next(r for r in identity if r['release']==release and r['construct']=='ALL')
        assert abs(current['repeated_site_mean_Brier']-float(siteold['Brier_observed']))<TOL
    outputs={'FILTER_CLASS_METRICS.tsv':classrows,'FILTER_CURVES.tsv':curve,'FIXED_CONSTRUCT_WEIGHTS.tsv':fixed,'LOSS_DECOMPOSITION.tsv':decompositions,'SITE_VARIANCE_IDENTITY.tsv':identity,'CONTROLLED_HIDING_BOUNDS.tsv':boundrows,'EXCLUSION_REASONS.tsv':reasonrows,'PARENT_RECONCILIATION.tsv':reconciliation}
    for name,rows in outputs.items(): write_tsv(out/name,rows)
    write_json(out/'RUN.json',dict(issue=protocol['issue'],created_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-started,cpu_seconds=time.process_time()-cpu,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,python=platform.python_version(),protocol_sha256=freeze['protocol_sha256'],code_sha256=sha(Path(__file__)),inputs=counts,output_rows={k:len(v) for k,v in outputs.items()},status='COMPLETED_TARGETED_CHECKS_PASS',boundary='Posthoc descriptive results; not new biological replication, native omission identification or journal readiness'))

if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('--workspace',required=True,type=Path); a.add_argument('--output',required=True,type=Path); x=a.parse_args(); run(x.workspace.resolve(),x.output)
