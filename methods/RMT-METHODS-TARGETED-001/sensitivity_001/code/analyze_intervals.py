"""ML interval sensitivity on frozen saved joins; no BAM or network access."""
import argparse,collections,csv,functools,gzip,hashlib,json,math,platform,resource,time
from datetime import datetime,timezone
from fractions import Fraction
from pathlib import Path

TOL=1e-10
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,x):
    with p.open('x') as f:json.dump(x,f,indent=2,allow_nan=False)
def table(p,rows):
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with p.open('x',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def mean(xs):return sum(xs)/len(xs) if xs and all(v is not None for v in xs) else None
def iprod(a,b):
    v=[x*y for x in a for y in b];return min(v),max(v)

@functools.lru_cache(None)
def bin_geometry(m,t):
    """Exact rational membership; returned loss endpoints are inf/sup envelopes."""
    lo,hi=Fraction(m,256),Fraction(m+1,256); q=Fraction(2*m+1,512);t=Fraction(str(t));b=1-t
    pieces=[]
    right=min(hi,b)
    if lo<right or (lo==right and lo<hi and lo<=b):pieces.append((lo,right))
    left=max(lo,t)
    if left<hi:pieces.append((left,hi))
    mandatory=t==Fraction(1,2) or lo>=t or hi<=b
    midpoint=max(q,1-q)>=t
    return float(lo),float(hi),midpoint,mandatory,pieces

def loss_extrema(pieces,y):
    vals=[float((x-y)**2) for piece in pieces for x in piece]
    return min(vals),max(vals)

def ratio_extreme(mandatory_n,mandatory_loss,optional,maximize=False):
    """Exact discrete count optimum apart from floating point arithmetic.

    Within one constant-loss optional bin, the linear-fractional objective is
    monotone (or constant) in k. Sorted endpoint bins admit a prefix optimum.
    An empty mandatory set returns extrema conditional on a nonempty selection.
    """
    candidates=sorted([(loss,n) for loss,n in optional if n],reverse=maximize)
    if not mandatory_n:
        return candidates[0][0] if candidates else None
    n=mandatory_n;s=mandatory_loss
    for loss,count in candidates:
        improves=(loss>s/n) if maximize else (loss<s/n)
        if improves:s+=loss*count;n+=count
        else:break
    return s/n

def envelopes(hist,t):
    fixed_n=mandatory_n=possible_n=0;midloss=flo=fhi=mlo=mhi=0.;optional_lo=[];optional_hi=[]
    for (y,m),n in sorted(hist.items()):
        lo,hi,mid,mandatory,pieces=bin_geometry(m,t)
        if mid:
            fixed_n+=n;midloss+=n*((m+.5)/256-y)**2
            a,b=loss_extrema([(Fraction(m,256),Fraction(m+1,256))],y);flo+=n*a;fhi+=n*b
        if not pieces:continue
        possible_n+=n;a,b=loss_extrema(pieces,y)
        if mandatory:mandatory_n+=n;mlo+=n*a;mhi+=n*b
        else:optional_lo.append((a,n));optional_hi.append((b,n))
    r={'midpoint_S':fixed_n,'minimum_S':mandatory_n,'maximum_S':possible_n,'ambiguous_S':possible_n-mandatory_n,'midpoint_Brier':midloss/fixed_n if fixed_n else None,'fixed_inf':flo/fixed_n if fixed_n else None,'fixed_sup':fhi/fixed_n if fixed_n else None,'variable_inf':ratio_extreme(mandatory_n,mlo,optional_lo),'variable_sup':ratio_extreme(mandatory_n,mhi,optional_hi,True),'fixed_guaranteed_nonempty':bool(fixed_n),'variable_guaranteed_nonempty':bool(mandatory_n),'variable_possible_nonempty':bool(possible_n)}
    assert mandatory_n<=fixed_n<=possible_n
    for layer in ('fixed','variable'):
        if r['midpoint_Brier'] is not None:assert r[layer+'_inf']-TOL<=r['midpoint_Brier']<=r[layer+'_sup']+TOL
    return r

def direction(a,b,guaranteed):
    if not guaranteed:return 'POTENTIALLY_UNDEFINED_REQUIRED_STRATUM',None,None
    lo,hi=b[0]-a[1],b[1]-a[0]
    return ('CERTIFIED_INCREASE' if lo>TOL else 'CERTIFIED_DECREASE' if hi<-TOL else 'NOT_CERTIFIED_BY_MARGINAL_ENVELOPES'),lo,hi

def decomposition_bounds(base,chosen,layer,t):
    if not all(r[layer+'_guaranteed_nonempty'] for r in base+chosen):return {'state':'POTENTIALLY_UNDEFINED_REQUIRED_CLASS'}
    if float(t)==.5:return dict(state='BASELINE_IDENTITY',within_inf=0.,within_sup=0.,weight_inf=0.,weight_sup=0.,sign_within='ZERO_IDENTITY',sign_weight='ZERO_IDENTITY',dominance='ZERO_IDENTITY')
    counts=[(r['midpoint_S'],r['midpoint_S']) if layer=='fixed' else (r['minimum_S'],r['maximum_S']) for r in chosen]
    n0=sum(r['midpoint_S'] for r in base);within=[0.,0.];weight=[0.,0.]
    for y in (0,1):
        a,b=counts[y],counts[1-y];w0=base[y]['midpoint_S']/n0
        w1=(a[0]/(a[0]+b[1]),a[1]/(a[1]+b[0]))
        l0=(base[y][layer+'_inf'],base[y][layer+'_sup']);l1=(chosen[y][layer+'_inf'],chosen[y][layer+'_sup'])
        w=iprod(((w1[0]+w0)/2,(w1[1]+w0)/2),(l1[0]-l0[1],l1[1]-l0[0]))
        z=iprod(((l1[0]+l0[0])/2,(l1[1]+l0[1])/2),(w1[0]-w0,w1[1]-w0))
        for j in (0,1):within[j]+=w[j];weight[j]+=z[j]
    def sign(v):return 'NEGATIVE_CERTIFIED' if v[1]<-TOL else 'POSITIVE_CERTIFIED' if v[0]>TOL else 'NOT_CERTIFIED'
    def absbounds(v):return (0. if v[0]<=0<=v[1] else min(abs(x) for x in v),max(abs(x) for x in v))
    aw,az=absbounds(within),absbounds(weight)
    return dict(state='CONSERVATIVE_COMPONENT_ENCLOSURES',within_inf=within[0],within_sup=within[1],weight_inf=weight[0],weight_sup=weight[1],sign_within=sign(within),sign_weight=sign(weight),dominance='WITHIN_ABSOLUTE_MAGNITUDE_CERTIFIED' if aw[0]>az[1]+TOL else 'WEIGHT_ABSOLUTE_MAGNITUDE_CERTIFIED' if az[0]>aw[1]+TOL else 'NOT_CERTIFIED')

def run(workspace,out):
    root=Path(__file__).resolve().parents[1];parent=root.parent
    p=json.loads((root/'PROTOCOL.json').read_text());pre=json.loads((root/'INPUTS.json').read_text())
    assert sha(root/'PROTOCOL.json')==pre['protocol_sha256']
    assert sha(parent/'run_001/results/FILTER_CURVES.tsv')==pre['source_curves_sha256']
    assert sha(parent/'run_001/results/LOSS_DECOMPOSITION.tsv')==pre['source_decomposition_sha256']
    out.mkdir(parents=True,exist_ok=False);start=time.monotonic();cpu=time.process_time()
    hists=collections.defaultdict(collections.Counter);counts=collections.defaultdict(collections.Counter)
    for source in pre['record_inputs']:
        path=workspace/source['path'];assert sha(path)==source['sha256'];release='2025_DRACH' if source['role'].startswith('development') else '2026_all5mer';n=0
        with gzip.open(path,'rt',newline='') as inp:
            for r in csv.DictReader(inp,delimiter='\t'):
                n+=1;y=int(r['truth_label']);assert y in (0,1);key=(release,r['contig'],y);counts[key]['opportunities']+=1
                eligible=r['alignment_status']=='MATCHED_READ_SITE';counts[key]['N']+=eligible
                if not r['p_mid']:continue
                m=int(r['ml_integer']);assert eligible and r['score_status']=='explicit' and 0<=m<=255
                assert r['score_kind']=='read_site_ML_interval' and r['modification_code'] in ('a','28871')
                assert float(r['p_mid'])==(m+.5)/256 and float(r['p_low'])==m/256 and float(r['p_high'])==(m+1)/256
                hists[key][(y,m)]+=1;counts[key]['O']+=1
        print(source['role'],n,flush=True)
    prior_file=parent/'run_001/results/FILTER_CLASS_METRICS.tsv';assert sha(prior_file)==pre['source_classes_sha256']
    with prior_file.open() as inp:prior={(r['release'],r['construct'],int(r['truth_class']),r['threshold']):r for r in csv.DictReader(inp,delimiter='\t') if r['policy']=='symmetric_confidence'}
    classes=[];summaries=[];fixed=[];decomps=[];recon=[];directions=[];lookup={};summary_lookup={};fixed_lookup={}
    releases=sorted({k[0] for k in hists});thresholds=p['thresholds']
    for release in releases:
        constructs=sorted({k[1] for k in hists if k[0]==release})
        for y in (0,1):
            for c in constructs:hists[(release,'ALL',y)].update(hists[(release,c,y)]);counts[(release,'ALL',y)].update(counts[(release,c,y)])
        for t in thresholds:
            for c in constructs+['ALL']:
                pair=[];both=collections.Counter()
                for y in (0,1):
                    k=(release,c,y);r=envelopes(hists[k],t);lookup[(release,c,y,t)]=r;pair.append(r);both.update(hists[k])
                    record=dict(release=release,construct=c,truth_class=y,threshold=t,analysis_class=p['analysis_class'],**dict(counts[k]),**r);classes.append(record)
                    old=prior[(release,c,y,str(t))]
                    for newkey,oldkey in [('midpoint_Brier','brier'),('midpoint_S','S'),('minimum_S','definite_retained'),('maximum_S','possible_retained')]:
                        a=r[newkey];b=float(old[oldkey]) if old[oldkey] else None;match=(a is None and b is None) or (a is not None and b is not None and abs(a-b)<=TOL)
                        assert match,(release,c,y,t,newkey,a,b);recon.append(dict(release=release,construct=c,truth_class=y,threshold=t,metric=newkey,old=b,new=a,match=match))
                obs=envelopes(both,t);s=dict(release=release,construct=c,threshold=t,analysis_class=p['analysis_class'],N=sum(counts[(release,c,y)]['N'] for y in (0,1)),O=sum(counts[(release,c,y)]['O'] for y in (0,1)),midpoint_S=obs['midpoint_S'],minimum_S=obs['minimum_S'],maximum_S=obs['maximum_S'],midpoint_observed=obs['midpoint_Brier'],midpoint_equal_class=mean([r['midpoint_Brier'] for r in pair]))
                for layer in ('fixed','variable'):
                    for end in ('inf','sup'):s[f'{layer}_observed_{end}']=obs[f'{layer}_{end}'];s[f'{layer}_equal_class_{end}']=mean([r[f'{layer}_{end}'] for r in pair])
                    s[f'{layer}_observed_guaranteed']=obs[layer+'_guaranteed_nonempty'];s[f'{layer}_equal_class_guaranteed']=all(r[layer+'_guaranteed_nonempty'] for r in pair)
                summaries.append(s);summary_lookup[(release,c,t)]=s
            ff=dict(release=release,threshold=t,K=len(constructs),analysis_class=p['analysis_class'],midpoint_equal_construct_class=mean([lookup[(release,c,y,t)]['midpoint_Brier'] for c in constructs for y in (0,1)]))
            for layer in ('fixed','variable'):
                cells=[lookup[(release,c,y,t)] for c in constructs for y in (0,1)]
                for end in ('inf','sup'):ff[f'{layer}_{end}']=mean([r[f'{layer}_{end}'] for r in cells])
                missing=[f'{c}:{y}' for c in constructs for y in (0,1) if not lookup[(release,c,y,t)][layer+'_guaranteed_nonempty']]
                ff[layer+'_guaranteed']=not missing;ff[layer+'_potentially_empty_cells']='|'.join(missing)
            fixed.append(ff);fixed_lookup[(release,t)]=ff
            for layer in ('fixed','variable'):
                d=decomposition_bounds([lookup[(release,'ALL',y,.5)] for y in (0,1)],[lookup[(release,'ALL',y,t)] for y in (0,1)],layer,t)
                decomps.append(dict(release=release,threshold=t,layer=layer,**d))
        pairs=sorted(set(zip(thresholds[:-1],thresholds[1:]))|{(.5,t) for t in thresholds[1:]})
        for a,b in pairs:
            for c in constructs+['ALL']:
                x,z=summary_lookup[(release,c,a)],summary_lookup[(release,c,b)]
                for layer in ('fixed','variable'):
                    for metric in ('observed','equal_class'):
                        status,lo,hi=direction((x[f'{layer}_{metric}_inf'],x[f'{layer}_{metric}_sup']),(z[f'{layer}_{metric}_inf'],z[f'{layer}_{metric}_sup']),x[f'{layer}_{metric}_guaranteed'] and z[f'{layer}_{metric}_guaranteed'])
                        directions.append(dict(release=release,construct=c,layer=layer,weighting=metric,threshold_from=a,threshold_to=b,delta_conservative_inf=lo,delta_conservative_sup=hi,status=status,interpretation='Sufficient marginal-envelope certificate; non-certificate does not establish reversal'))
            for layer in ('fixed','variable'):
                x,z=fixed_lookup[(release,a)],fixed_lookup[(release,b)];status,lo,hi=direction((x[layer+'_inf'],x[layer+'_sup']),(z[layer+'_inf'],z[layer+'_sup']),x[layer+'_guaranteed'] and z[layer+'_guaranteed'])
                directions.append(dict(release=release,construct='FIXED_ALL_CONSTRUCTS',layer=layer,weighting='equal_construct_class',threshold_from=a,threshold_to=b,delta_conservative_inf=lo,delta_conservative_sup=hi,status=status,interpretation='Fixed original K; sufficient certificate only'))
    # Bind midpoint decomposition checks to the unchanged prior class decomposition.
    with (parent/'run_001/results/LOSS_DECOMPOSITION.tsv').open() as inp:olddecomp={(r['release'],float(r['threshold'])):r for r in csv.DictReader(inp,delimiter='\t') if r['stratum_scheme']=='truth_class' and r['policy']=='symmetric_confidence'}
    for d in decomps:
        old=olddecomp[(d['release'],d['threshold'])]
        for field,oldfield in [('within','within_loss_component'),('weight','weight_component')]:
            point=float(old[oldfield]);d['midpoint_'+field]=point
            if field+'_inf' in d:assert d[field+'_inf']-TOL<=point<=d[field+'_sup']+TOL
    outputs={'CLASS_INTERVALS.tsv':classes,'SUMMARY_INTERVALS.tsv':summaries,'FIXED_CONSTRUCT_INTERVALS.tsv':fixed,'DIRECTION_CERTIFICATION.tsv':directions,'DECOMPOSITION_ENCLOSURES.tsv':decomps,'PARENT_RECONCILIATION.tsv':recon}
    for name,rows in outputs.items():table(out/name,rows)
    dump(out/'RUN.json',dict(status='COMPLETED_WITH_TARGETED_CHECKS_PASS',created_utc=datetime.now(timezone.utc).isoformat(),wall_seconds=time.monotonic()-start,cpu_seconds=time.process_time()-cpu,max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,threads=1,python=platform.python_version(),protocol_sha256=sha(root/'PROTOCOL.json'),code_sha256=sha(Path(__file__)),row_counts={n:len(r) for n,r in outputs.items()},scope='Posthoc marginal ML-model-score sensitivity only; not native omission identification or biological truth'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--workspace',required=True,type=Path);p.add_argument('--output',required=True,type=Path);a=p.parse_args();run(a.workspace,a.output)
