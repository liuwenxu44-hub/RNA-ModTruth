"""Streaming posthoc evaluation of saved joins. No BAM, networking or caller I/O."""
import argparse, collections, csv, gzip, hashlib, heapq, json, math, resource, time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).resolve().parents[1]
E=ROOT/'execution/RMT-EXEC-001'
P=json.loads((HERE/'PROTOCOL.json').read_text())
SEEDS=P['sampling']['seeds']
CONDS=['canonical','site_mean_as_read','missing_zero','high_p_filter']

def sha(p):
    with open(p,'rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()
def dump(p,obj):
    with p.open('x') as f: json.dump(obj,f,indent=2,allow_nan=False)
def tsv(p,rows):
    with p.open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
def records(p):
    with gzip.open(p,'rt',newline='') as f: yield from csv.DictReader(f,delimiter='\t')
def rid(r): return (r['sample_id'],r['run_id'],r['contig'],r['read_id'])
def safe(x): return x if math.isfinite(x) else 'INFINITE'
def loss(p,y):
    q=p if y else 1-p
    return -math.log(q) if q>0 else math.inf

class Site:
    def __init__(self): self.opp=0; self.eligible=0; self.hist=collections.Counter(); self.missing=collections.Counter()

class Group:
    def __init__(self):
        self.den=collections.Counter(); self.n=[0,0]; self.b=[0.,0.]; self.ll=[0.,0.]
        self.blo=[0.,0.]; self.bhi=[0.,0.]; self.llo=[0.,0.]; self.llhi=[0.,0.]
        self.invalid=[0,0]; self.conf=collections.Counter(TP=0,TN=0,FP=0,FN=0)
    def add(self,y,p,n,lo=None,hi=None):
        if not n:return
        self.n[y]+=n; self.b[y]+=n*(p-y)**2; self.ll[y]+=n*loss(p,y)
        self.conf['TP' if y and p>=.5 else 'FN' if y else 'FP' if p>=.5 else 'TN']+=n
        if lo is None: self.invalid[y]+=n
        else:
            self.blo[y]+=n*min((lo-y)**2,(hi-y)**2); self.bhi[y]+=n*max((lo-y)**2,(hi-y)**2)
            self.llo[y]+=n*min(loss(lo,y),loss(hi,y)); self.llhi[y]+=n*max(loss(lo,y),loss(hi,y))
    def row(self,release,subset,condition,construct):
        n=sum(self.n); b=[self.b[y]/self.n[y] if self.n[y] else None for y in (0,1)]
        allnative=n and not sum(self.invalid)
        tp,tn,fp,fn=(self.conf[k] for k in ('TP','TN','FP','FN'))
        tpr=tp/(tp+fn) if tp+fn else None; fpr=fp/(fp+tn) if fp+tn else None
        d=dict(release=release,subset=subset,condition=condition,construct=construct,
          analysis_class='POSTHOC_ROBUSTNESS_ANALYSIS',evidence_class='OBSERVED_NATIVE_INPUT' if condition=='canonical' else 'CONTROLLED_ERROR_INJECTION',
          truth_unit='provider_construct_condition',metric_unit='read_site',
          opportunities=self.den['opp'],eligible=self.den['eligible'],native_explicit=self.den['explicit'],retained_records=self.den['retained'],
          opportunities_control=self.den['opp_0'],opportunities_modified=self.den['opp_1'],
          eligible_control=self.den['eligible_0'],eligible_modified=self.den['eligible_1'],
          explicit_control=self.den['explicit_0'],explicit_modified=self.den['explicit_1'],
          retained_control=self.den['retained_0'],retained_modified=self.den['retained_1'],
          scored=n,scored_control=self.n[0],scored_modified=self.n[1],
          Brier_control=b[0],Brier_modified=b[1],Brier_observed=sum(self.b)/n if n else None,
          Brier_equal_class=sum(b)/2 if all(v is not None for v in b) else None,
          log_loss=safe(sum(self.ll)/n) if n else None,
          log_loss_control=safe(self.ll[0]/self.n[0]) if self.n[0] else None,
          log_loss_modified=safe(self.ll[1]/self.n[1]) if self.n[1] else None,
          balanced_accuracy=(tpr+1-fpr)/2 if tpr is not None and fpr is not None else None,FPR=fpr,
          **dict(self.conf),native_interval_points=n-sum(self.invalid),
          Brier_interval_infimum=sum(self.blo)/n if allnative else None,
          Brier_interval_supremum=sum(self.bhi)/n if allnative else None,
          log_loss_interval_infimum=safe(sum(self.llo)/n) if allnative else None,
          log_loss_interval_supremum=safe(sum(self.llhi)/n) if allnative else None,
          conditional_scope='eligible explicit probability subset' if condition=='canonical' else 'invalid unit: repeated site means' if condition=='site_mean_as_read' else 'invalid probability: missing replaced by zero' if condition=='missing_zero' else 'selected p>=0.9 subset; whole-population claim invalid')
        for y in (0,1):
            d['explicit_coverage_'+str(y)]=self.den['explicit_'+str(y)]/self.den['opp_'+str(y)] if self.den['opp_'+str(y)] else None
        return d

def run(out):
    t=time.monotonic(); cpu=time.process_time(); out.mkdir(parents=True,exist_ok=False)
    frozen=json.loads((HERE/'INPUTS.json').read_text())
    assert sha(HERE/'PROTOCOL.json')==frozen['protocol_sha256']
    rows_out=[]; selected_output=[]; frame_output=[]; counts_output=[]; reconciliation=[]
    allsites={}
    for inp in frozen['record_inputs']:
        p=ROOT/inp['path']; assert sha(p)==inp['sha256']
        release='2025_DRACH' if inp['role'].startswith('development') else '2026_all5mer'
        seen=collections.defaultdict(set); old_reads=collections.defaultdict(set); fullcount=collections.Counter()
        # First pass uses identifiers and inclusion flags only; no probability in selection.
        for r in records(p):
            key=rid(r); seen[key[:3]].add(key)
            if r['selected_for_variants']=='1': old_reads[key[:3]].add(key)
            fullcount['opportunities']+=1
            fullcount['eligible']+=r['alignment_status']=='MATCHED_READ_SITE'
            fullcount['explicit']+=bool(r['p_mid'])
            fullcount['selected_original_opportunities']+=r['selected_for_variants']=='1'
            fullcount['status_'+r['score_status']]+=1
        chosen=collections.defaultdict(list)
        for group,keys in sorted(seen.items()):
            frame_output.append(dict(release=release,sample_id=group[0],run_id=group[1],construct=group[2],available_overlap_reads=len(keys),original_selected_reads=len(old_reads[group])))
            for seed in SEEDS:
                def rank(k): return (hashlib.sha256(json.dumps([seed,*k],separators=(',',':'),ensure_ascii=False).encode()).hexdigest(),k)
                chosen_keys=heapq.nsmallest(100,keys,key=rank)
                for key in chosen_keys:
                    chosen[key].append('seed_'+str(seed))
                    selected_output.append(dict(release=release,seed=seed,sample_id=key[0],run_id=key[1],construct=key[2],read_id=key[3],rank_sha256=rank(key)[0]))
        del seen,old_reads
        # Compress to exact ML histograms within sample/site/subset; no raw-row copies.
        for r in records(p):
            subsets=['full']+chosen.get(rid(r),[])
            if r['selected_for_variants']=='1': subsets.append('original_first100')
            eligible=r['alignment_status']=='MATCHED_READ_SITE'; explicit=bool(r['p_mid'])
            if explicit:
                ml=int(r['ml_integer'])
                assert eligible and float(r['p_mid'])==(ml+.5)/256 and r['score_status']=='explicit'
                assert float(r['p_low'])==ml/256 and float(r['p_high'])==(ml+1)/256
            y=int(r['truth_label']); assert y in (0,1)
            for subset in subsets:
                k=(release,subset,r['sample_id'],r['contig'],y,r['ref_pos0'],r['strand'])
                s=allsites.setdefault(k,Site()); s.opp+=1; s.eligible+=eligible
                if explicit: s.hist[ml]+=1
                else: s.missing[r['score_status']]+=1
        summary=json.loads(p.with_name('summary.json').read_text())['counts']
        expected=dict(opportunities=summary['truth_overlapping_alignment_sites'],explicit=summary['evaluable_explicit_read_sites'],selected_original_opportunities=summary['selected_variant_sites'])
        assert all(fullcount[k]==v for k,v in expected.items()), (inp['role'],dict(fullcount),expected)
        counts_output.append(dict(input_role=inp['role'],source=inp['path'],**dict(fullcount)))
        print(inp['role'],dict(fullcount),flush=True)
    groups=collections.defaultdict(Group)
    for (release,subset,sample,construct,y,pos,strand),s in allsites.items():
        explicit=sum(s.hist.values()); filt=sum(n for m,n in s.hist.items() if (m+.5)/256>=.9)
        for cond in CONDS:
            for unit in (construct,'ALL'):
                g=groups[(release,subset,cond,unit)]
                for k,v in [('opp',s.opp),('eligible',s.eligible),('explicit',explicit),('retained',filt if cond=='high_p_filter' else s.opp)]:
                    g.den[k]+=v; g.den[k+'_'+str(y)]+=v
                if cond=='site_mean_as_read':
                    if explicit: g.add(y,sum((m+.5)/256*n for m,n in s.hist.items())/explicit,explicit)
                else:
                    for ml,n in s.hist.items():
                        p=(ml+.5)/256
                        if cond=='high_p_filter' and p<.9: continue
                        g.add(y,p,n,ml/256,(ml+1)/256)
                    if cond=='missing_zero': g.add(y,0.,s.eligible-explicit)
    for key,g in sorted(groups.items()): rows_out.append(g.row(*key))
    # Independent legacy neutral-consumer values must match for original subsets.
    for release,resultroot in [('2025_DRACH',E/'results/development_test_final'),('2026_all5mer',E/'results/evaluation_run_v1/results/evaluation')]:
        for cond,case in zip(CONDS,['00_canonical','07_site_score_as_read_probability','09_missing_probability_imputed_zero','08_filtered_output_claimed_all_records']):
            original=json.loads((resultroot/'metrics'/f'{case}.json').read_text())
            r=next(r for r in rows_out if (r['release'],r['subset'],r['condition'],r['construct'])==(release,'original_first100',cond,'ALL'))
            pairs=[('Brier_observed',original['Brier']),('log_loss',original['log_loss']),('balanced_accuracy',original['balanced_accuracy']),('FPR',original['FPR']),('scored',original['counts']['evaluable_points'])]
            okay=all(r[k]==v if not isinstance(v,float) else math.isclose(r[k],v,rel_tol=1e-12,abs_tol=1e-12) for k,v in pairs)
            reconciliation.append(dict(release=release,condition=cond,source=str((resultroot/'metrics'/f'{case}.json').relative_to(ROOT)),status='PASS' if okay else 'FAIL'))
    assert all(r['status']=='PASS' for r in reconciliation), reconciliation
    dump(out/'METRICS.json',rows_out); tsv(out/'METRICS.tsv',rows_out)
    tsv(out/'SELECTED_READS_LOCAL_ONLY.tsv',selected_output); tsv(out/'SAMPLING_FRAME.tsv',frame_output)
    dump(out/'NATIVE_COUNTS.json',counts_output); dump(out/'LEGACY_RECONCILIATION.json',reconciliation)
    dump(out/'RUN.json',dict(protocol_sha256=sha(HERE/'PROTOCOL.json'),code_sha256=sha(Path(__file__)),wall_seconds=time.monotonic()-t,cpu_seconds=time.process_time()-cpu,peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,metric_rows=len(rows_out),sampled_memberships=len(selected_output),status='PASS',analysis_class='POSTHOC_ROBUSTNESS_ANALYSIS',network=False,threads=1))
    print('COMPLETE',len(rows_out),'metric rows; all eight legacy comparisons passed',flush=True)

if __name__=='__main__':
    a=argparse.ArgumentParser(); a.add_argument('--output',required=True); args=a.parse_args(); run(Path(args.output))
