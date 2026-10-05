"""Formatting frozen outputs only; no new inference or selection."""
import csv, hashlib, json
from pathlib import Path

r=Path(__file__).resolve().parents[1]; p=r/'run_001/results'
def rows(name):
    with (p/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
def f(x):return 'NA' if x=='' else f'{float(x):.6f}'
lines=['# Targeted supplementary numeric results','', 'These are posthoc descriptive midpoint-based results on saved eligible explicit scores, not population-wide performance or a model ranking. All seven frozen thresholds are reported. Percentages refer to the named denominator; no threshold is selected as optimal.','', '## Complete symmetric-filter sequence','', '| Release | t | N eligible | O explicit | S retained | S/O | S/N | Observed Brier | Equal-class Brier | Equal-class risk | Ambiguous ML retention |','|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for x in rows('FILTER_CURVES.tsv'):
    if x['construct']=='ALL' and x['policy']=='symmetric_confidence':
        v=[x['release'],x['threshold'],x['N'],x['O'],x['S'],f(x['coverage_S_over_O']),f(x['coverage_S_over_N']),f(x['Brier_observed']),f(x['Brier_equal_class']),f(x['risk_equal_class']),x['ambiguous_retention']]
        lines.append('| '+' | '.join(v)+' |')
lines+=['','## Class decomposition: prespecified reference and illustrative high thresholds','', 'Complete class and construct decompositions are in LOSS_DECOMPOSITION.tsv; these rows illustrate interpretation, not cutoff selection. A positive within-class component means increased loss in this symmetric algebra, not causal deterioration.','', '| Release | Policy | t | Selected − baseline Brier | Within-class component | Class-weight component | Residual |','|---|---|---:|---:|---:|---:|---:|']
for x in rows('LOSS_DECOMPOSITION.tsv'):
    if x['stratum_scheme']=='truth_class' and x['threshold'] in ('0.9','0.99'):
        lines.append('| '+' | '.join([x['release'],x['policy'],x['threshold']]+[f(x[k]) for k in ['delta','within_loss_component','weight_component']]+[x['addback_residual']])+' |')
lines+=['','## Uniform-label site identity','', '| Release | Explicit read-sites | Sample-site groups | Mixed-label groups | Read Brier | Repeated-mean Brier | Within-site population variance |','|---|---:|---:|---:|---:|---:|---:|']
for x in rows('SITE_VARIANCE_IDENTITY.tsv'):
    if x['construct']=='ALL':lines.append('| '+' | '.join([x['release'],x['explicit_records'],x['site_groups'],x['mixed_label_groups']]+[f(x[k]) for k in ['read_Brier','repeated_site_mean_Brier','population_within_site_variance']])+' |')
lines+=['','The 40 and 512 groups are sample/run/construct/site/strand groups, not independent biological replicates. The repeated site mean is not recovered read probability.','', '## Controlled hiding only','', '| Release | Weighting | Original O | Hidden | Loss infimum | Known midpoint loss | Loss supremum |','|---|---|---:|---:|---:|---:|---:|']
for x in rows('CONTROLLED_HIDING_BOUNDS.tsv'):
    if x['truth_class'] in ('ALL','EQUAL_CLASS'):lines.append('| '+' | '.join([x['release'],x['truth_class'],x['N'],x['hidden']]+[f(x[k]) for k in ['loss_infimum','known_midpoint_loss','loss_supremum']])+' |')
lines+=['','Scores were hidden only when their original ML upper endpoint was ≤0.25. These envelopes are a controlled implementation check, not identified native omission bounds, confidence intervals, or occupancy estimates.','']
with (r/'NUMERIC_RESULTS.md').open('x') as out:out.write('\n'.join(lines))
print('Created NUMERIC_RESULTS.md from frozen TSVs')
