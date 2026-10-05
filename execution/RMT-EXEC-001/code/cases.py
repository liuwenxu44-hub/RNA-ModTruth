"""Frozen, outcome-independent controlled transformations; never imports a checker."""
import argparse
import collections
import copy
import json
from pathlib import Path
from bundle import load, digest, write_json, write_table, write_fasta

ERRORS = ['sample_reference_BED_mismatch','declared_zero_one_based_coordinate_error',
          'strand_or_reference_version_conflict','lost_sample_join_key',
          'read_index_reuse_across_runs','duplicate_read_site_counting',
          'site_score_as_read_probability','filtered_output_claimed_all_records',
          'missing_probability_imputed_zero','modification_code_confusion']
LEGAL = ['row_order','compression','irrelevant_column','documented_bijective_rename',
         'disjoint_merge','documented_coordinate_conversion']


def generate(parent, output):
    parent, output = Path(parent), Path(output)
    output.mkdir(parents=True, exist_ok=False)
    pc, pr, pt, pf, pm = load(parent)
    for row in pr:
        row.pop('_locator')
    parents = {n:digest(parent/n) for n in ('contract.json','predictions.tsv.gz','truth.tsv','reference.fa','read_map.tsv')}
    records = []
    for index, family in enumerate(['canonical']+ERRORS+LEGAL):
        contract, rows, truth, refs, mapping = copy.deepcopy((pc,pr,pt,pf,pm))
        path = output/f'{index:02d}_{family}'
        path.mkdir()
        affected, changes = [], []
        if family == 'sample_reference_BED_mismatch':
            samples = [s['sample_id'] for s in contract['samples']]
            swap = dict(zip(samples, reversed(samples)))
            for row in truth:
                row['sample_id'] = swap[row['sample_id']]
            affected = [r['evidence_id'] for r in rows]
            changes = ['Swap control/modified truth sample assignment; preserve original BED labels, source SHA and line.']
        elif family in ('declared_zero_one_based_coordinate_error','documented_coordinate_conversion'):
            for row in rows:
                row['position'] = str(int(row['position'])+1)
            contract['producer_origin'] = 1
            if family == 'documented_coordinate_conversion':
                contract['consumer_origin'] = 1
            affected = [r['evidence_id'] for r in rows]
            changes = ['Prediction coordinates +1; producer origin 1; consumer '+str(contract['consumer_origin'])]
        elif family == 'strand_or_reference_version_conflict':
            rows[0]['strand'] = '-' if rows[0]['strand'] == '+' else '+'
            affected = [rows[0]['evidence_id']]
            changes = ['Invert first observation strand, not source alignment or BED.']
        elif family == 'lost_sample_join_key':
            contract['truth_join_keys'].remove('sample_id')
            affected = [r['evidence_id'] for r in rows]
            changes = ['Drop sample_id from downstream truth join keys; rows still retain verifiable identities.']
        elif family == 'read_index_reuse_across_runs':
            target = rows[0]
            other = next(r for r in mapping if r['read_index'] == target['read_index'] and r['sample_id'] != target['sample_id'])
            target.update({k:other[k] for k in ('sample_id','run_id','read_id')})
            affected = [target['evidence_id']]
            changes = ['Resolve one local adapter read index using other sample/run namespace; preserve BAM ordinal/evidence ID.']
        elif family == 'duplicate_read_site_counting':
            target = next(r for r in rows if r['score_status'] == 'explicit')
            rows.append(copy.deepcopy(target))
            affected = [target['evidence_id']]
            changes = ['Append one fully identical source-linked explicit-probability observation.']
        elif family == 'site_score_as_read_probability':
            by_site = collections.defaultdict(list)
            for row in rows:
                if row['p_mid']:
                    by_site[(row['sample_id'],row['contig'],row['position'])].append(float(row['p_mid']))
            for row in rows:
                if row['p_mid']:
                    values = by_site[(row['sample_id'],row['contig'],row['position'])]
                    row.update(p_mid=str(sum(values)/len(values)),score_kind='site_aggregate_mean')
                    affected.append(row['evidence_id'])
            changes = ['Controlled aggregation of stored fixed probabilities to a site mean, then misuse as read probability; not m6Anet probabilities or biological truth.']
        elif family == 'filtered_output_claimed_all_records':
            rows = [r for r in rows if r['p_mid'] and float(r['p_mid']) >= .9]
            kept = {r['evidence_id'] for r in rows}
            affected = [r['evidence_id'] for r in pr if r['evidence_id'] not in kept]
            changes = ['Keep explicit probabilities >=0.9 without narrowing the selected-opportunity population claim.']
        elif family == 'missing_probability_imputed_zero':
            for row in rows:
                if not row['p_mid']:
                    row['p_mid'] = '0'
                    affected.append(row['evidence_id'])
            changes = ['Fill all missing probability points with zero; original ML and skip/missing states retained.']
        elif family == 'modification_code_confusion':
            for row in rows:
                row['modification_code'] = '17596'
            affected = [r['evidence_id'] for r in rows]
            changes = ['Relabel m6A as inosine without changing stored probabilities.']
        elif family == 'row_order':
            rows.reverse()
            changes = ['Reverse row order; unchanged logical multiset.']
        elif family == 'compression':
            contract['prediction_parts'] = ['predictions.tsv']
            changes = ['Lossless gzip-to-plain TSV.']
        elif family == 'irrelevant_column':
            for row in rows:
                row['export_note'] = 'ignored provenance-neutral column'
            changes = ['Add an irrelevant column without changing required fields.']
        elif family == 'documented_bijective_rename':
            rename = {k:f'alias_{i:03d}' for i,k in enumerate(refs)}
            contract['contig_aliases'] = {v:k for k,v in rename.items()}
            refs = {rename[k]:v for k,v in refs.items()}
            for row in rows+truth:
                row['contig'] = rename[row['contig']]
            changes = ['Bijective, explicitly documented reference alias mapping in predictions, truth and FASTA.']
        elif family == 'disjoint_merge':
            contract['prediction_parts'] = ['part_A.tsv.gz','part_B.tsv.gz']
            changes = ['Two disjoint row partitions consumed as a union; no repeated observation.']
        elif family == 'canonical':
            changes = ['Identity copy.']
        if len(contract['prediction_parts']) == 2:
            cut = len(rows)//2
            write_table(path/contract['prediction_parts'][0],rows[:cut])
            write_table(path/contract['prediction_parts'][1],rows[cut:])
        else:
            write_table(path/contract['prediction_parts'][0],rows)
        write_table(path/'truth.tsv',truth)
        write_table(path/'read_map.tsv',mapping)
        write_fasta(path/'reference.fa',refs)
        write_json(path/'contract.json',contract)
        entry = dict(case_id=path.name,family=family,role=contract['role'],
            kind='CONTROLLED_ERROR' if family in ERRORS else 'LEGAL_TRANSFORM' if family in LEGAL else 'CANONICAL',
            source='CONTROLLED_INJECTION_NOT_OBSERVED_DATABASE_ERROR',parent_files_sha256=parents,
            operations=changes,affected_evidence_ids=sorted(set(affected)),
            expected_constraint='REJECT_OR_SCOPE_UNSUPPORTED' if family in ERRORS else 'ACCEPT',
            generated_files_sha256={n:digest(path/n) for n in contract['prediction_parts']+['truth.tsv','reference.fa','read_map.tsv','contract.json']})
        # This truth about the injection is excluded from checker input; the evaluator reads it only afterwards.
        records.append(entry)
    write_json(output/'CASE_INDEX.json',records)
    print(json.dumps({'output':str(output),'cases':len(records),'role':pc['role']}))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--parent',required=True)
    p.add_argument('--output',required=True)
    a = p.parse_args()
    generate(a.parent,a.output)
