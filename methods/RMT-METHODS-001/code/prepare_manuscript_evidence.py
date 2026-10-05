"""Source-table extraction of completed results; no new experiment."""
import csv, hashlib, json, math, socket
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]; M=ROOT/'methods/RMT-METHODS-001'; E=ROOT/'execution/RMT-EXEC-001'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def table(p,rows):
 with p.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
 R=M/'results/server_run_001';out=M/'evidence/server_002';out.mkdir(parents=True,exist_ok=False)
 metrics=json.loads((R/'robustness/METRICS.json').read_text());facts=json.loads((R/'facts/FACTS.json').read_text());checks=json.loads((R/'facts/METHOD_CHECKS.json').read_text())
 coverage=[];composition=[];summary=[];contrasts=[];bounds=[]
 for r in metrics:
  if r['subset']=='full' and r['condition']=='canonical':
   for y,label in [(0,'Control'),(1,'Modified')]:
    suffix='control' if y==0 else 'modified';opp=r['opportunities_'+suffix];eligible=r['eligible_'+suffix];exp=r['explicit_'+suffix]
    coverage.append(dict(release=r['release'],construct=r['construct'],truth_class=label,opportunities=opp,eligible=eligible,explicit=exp,eligible_fraction=eligible/opp if opp else None,explicit_fraction=exp/opp if opp else None,explicit_given_eligible=exp/eligible if eligible else None,analysis_class=r['analysis_class']))
  if r['construct']=='ALL':
   summary.append(r)
   composition.append(dict(release=r['release'],subset=r['subset'],condition=r['condition'],opportunities=r['opportunities'],eligible=r['eligible'],native_explicit=r['native_explicit'],retained=r['retained_records'],scored=r['scored'],control=r['scored_control'],modified=r['scored_modified'],modified_fraction=r['scored_modified']/r['scored'] if r['scored'] else None,scored_opportunity_fraction=r['scored']/r['opportunities']))
   if r['condition']=='canonical':bounds.append({k:r[k] for k in ['release','subset','Brier_observed','Brier_interval_infimum','Brier_interval_supremum','log_loss','log_loss_interval_infimum','log_loss_interval_supremum','scored']})
  if r['condition']!='canonical':
   base=next(x for x in metrics if (x['release'],x['subset'],x['construct'],x['condition'])==(r['release'],r['subset'],r['construct'],'canonical'))
   for key in ['Brier_observed','Brier_equal_class','Brier_control','Brier_modified']:
    contrasts.append(dict(release=r['release'],subset=r['subset'],construct=r['construct'],condition=r['condition'],metric=key,canonical=base[key],transformed=r[key],delta=r[key]-base[key] if r[key] is not None and base[key] is not None else None,analysis_class='POSTHOC_ROBUSTNESS_ANALYSIS',evidence_class='CONTROLLED_ERROR_INJECTION'))
 table(out/'COVERAGE.tsv',coverage);table(out/'COMPOSITION.tsv',composition);table(out/'AGGREGATE_METRICS.tsv',summary);table(out/'CONSTRUCT_CONTRASTS.tsv',contrasts);table(out/'QUANTIZATION_BOUNDS.tsv',bounds)
 native=[]
 for release,rel in [('2025_DRACH','results/native_development_v1'),('2026_all5mer','results/evaluation_run_v1/results/native_evaluation')]:
  for mode in ['no_confidence_filter','default_confidence_filter']:
   p=E/rel/(mode+'.tsv');obj={k:json.loads(v) for k,v in (line.split(': ',1) for line in p.read_text().splitlines() if line)}
   for which in ['raw','filtered']:
    if which+'_contingency_table' not in obj:continue
    ct=obj[which+'_contingency_table'];n=sum(sum(row[1:]) for row in ct[1:])
    native.append(dict(release=release,mode=mode,output_population=which,metric='modkit_multiclass_accuracy_percent',accuracy=obj[which+'_accuracy'],denominator=n,control=sum(ct[1][1:]),modified=sum(ct[2][1:]),filter_threshold=obj['filter_threshold'],reported_percent_removed=obj.get('percent_of_mod_called_removed',0),scope='modkit class-balanced multiclass confusion table; not fixed binary explicit-m6A population',source=str(p.relative_to(ROOT)),source_sha256=sha(p),source_selector=which+'_accuracy;'+which+'_contingency_table',evidence_class='OBSERVED_NATIVE_INPUT',execution='HISTORICAL_EXEC_NOT_RERUN'))
 table(out/'NATIVE_MODKIT.tsv',native)
 methods=[];endpoint_source=[]
 for role in ['DEVELOPMENT','EVALUATION']:
  ci=E/('controlled_cases/development_v2/CASE_INDEX.json' if role=='DEVELOPMENT' else 'results/evaluation_run_v1/controlled_cases/evaluation/CASE_INDEX.json');kinds={x['case_id']:x['kind'] for x in json.loads(ci.read_text())}
  resultroot=E/('results/development_test_final' if role=='DEVELOPMENT' else 'results/evaluation_run_v1/results/evaluation')
  with (resultroot/'method_results.tsv').open() as f:endpoints=list(csv.DictReader(f,delimiter='\t'))
  endpoint_source.extend(endpoints)
  for method in ['B0','B1','B2','RMT']:
   rows=[r for r in checks if r['role']==role and r['method']==method];errors=[r for r in rows if kinds[r['case_id']]=='CONTROLLED_ERROR'];legal=[r for r in rows if kinds[r['case_id']]=='LEGAL_TRANSFORM']
   er=[r for r in endpoints if r['method']==method and r['kind']=='CONTROLLED_ERROR']
   unassessed=sum(r['endpoint_result']=='UNASSESSED' for r in er)
   assert unassessed==({'B0':10,'B1':3,'B2':0,'RMT':0}[method])
   methods.append(dict(role=role,method=method,controlled_error_families=len(errors),reported_detected=sum(r['reported_detects_error']=='True' for r in errors),scientific_unassessed=unassessed,legal_transforms=len(legal),legal_false_rejections=sum(r['reported_legal_false_rejection']=='True' for r in legal),safe_recoveries=sum(r.get('recovery',{}).get('written',False) for r in errors),status_scope='structural only' if method=='B0' else 'native tools and explicit joins; three export-semantic tasks unassessed' if method=='B1' else 'source-attested semantic checks',evidence_class='CONTROLLED_ERROR_INJECTION + LEGAL_TRANSFORMATION + SAFE_RECOVERY',endpoint_source=str((resultroot/'method_results.tsv').relative_to(ROOT)),endpoint_source_sha256=sha(resultroot/'method_results.tsv')))
 table(out/'METHOD_SCOPE.tsv',methods)
 table(out/'METHOD_ENDPOINT_SOURCE.tsv',endpoint_source)
 effects=[]
 for r in facts:
  if r['evidence_class']=='SAFE_RECOVERY':continue
  m=r['recalculated_metrics'];base=next(x['recalculated_metrics'] for x in facts if x['role']==r['role'] and x['case_id']=='00_canonical')
  effects.append(dict(role=r['role'],case_id=r['case_id'],evidence_class=r['evidence_class'],scored=m['counts']['evaluable_points'],Brier=m['Brier'],Brier_equal_class=r['Brier_equal_class'],balanced_accuracy=m['balanced_accuracy'],log_loss=m['log_loss'],delta_Brier=m['Brier']-base['Brier'] if m['Brier'] is not None else None,outcome='NO_RESULT' if m['Brier'] is None else 'NUMERICALLY_UNCHANGED_BRIER' if m['Brier']==base['Brier'] else 'CHANGED_BRIER',source_metric=r['source_metric'],source_metric_sha256=r['source_metric_sha256']))
 table(out/'ALL_CASE_EFFECTS.tsv',effects)
 notes=[]
 for rel in ['2025_DRACH','2026_all5mer']:
  full=next(r for r in summary if r['release']==rel and r['subset']=='full' and r['condition']=='canonical');orig=next(r for r in summary if r['release']==rel and r['subset']=='original_first100' and r['condition']=='canonical');seeds=[r['Brier_observed'] for r in summary if r['release']==rel and r['subset'].startswith('seed_') and r['condition']=='canonical']
  notes.append(dict(release=rel,original_Brier=orig['Brier_observed'],full_Brier=full['Brier_observed'],five_seed_min=min(seeds),five_seed_max=max(seeds),original_in_seed_range=min(seeds)<=orig['Brier_observed']<=max(seeds),scope='Five deterministic technical subsets, not a confidence interval'))
 with (out/'SUMMARY.json').open('x') as f:json.dump(dict(status='PASS',host=socket.gethostname(),analysis_input_sha256=sha(R/'robustness/METRICS.json'),new_statistical_tests=False,seed_summary=notes,undefined_equal_class_rows=sum(r['Brier_equal_class'] is None for r in metrics)),f,indent=2)
 print(json.dumps(notes,indent=2));print('SOURCE TABLE EXTRACTION COMPLETE',len(metrics),len(facts),len(checks))
if __name__=='__main__':main()
