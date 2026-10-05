"""RMT source-bound export validator. No error-family names or expected case labels.

The trusted ledger is re-verifiable from native files; it is not RMT decision gold.
This checker does not certify a novel biological label, model independence or purity.
"""
import argparse
import collections
import copy
from pathlib import Path
from bundle import load, read_json, digest, write_json, write_table, write_fasta

RECOVERABLE = {'OFFICIAL_TRUTH_MAPPING','DECLARED_COORDINATE_CONFLICT',
               'READ_IDENTITY_CONFLICT','EXACT_DUPLICATE'}


def read_maps_equal(observed, source):
    """v0.2 regression fix: orderless keys, no silent deduplication."""
    fields = ('sample_id', 'run_id', 'read_index')
    def indexed(rows):
        out = {}
        for row in rows:
            if any(k not in row or row[k] == '' for k in fields):
                return None
            key = tuple(row[k] for k in fields)
            if key in out:
                return None
            out[key] = dict(row)
        return out
    left, right = indexed(observed), indexed(source)
    return left is not None and right is not None and left == right


def validate(bundle, evidence):
    bundle, evidence = Path(bundle), Path(evidence)
    c, rows, truth, refs, mapping = load(bundle)
    ec, erows, etruth, erefs, emapping = load(evidence)
    links = read_json(evidence/'SOURCE_LINKS.json')
    if any(digest(evidence/name) != sha for name,sha in links['files'].items()):
        raise ValueError('TRUSTED_SOURCE_LEDGER_INTEGRITY_FAILURE')
    expected = {r['evidence_id']:r for r in erows}
    issues = []
    all_ids = sorted(expected)
    def issue(rule, ids, old, new, locator, evidence_locator):
        issues.append(dict(rule=rule, evidence_ids=sorted(set(ids)), observed=old,
            expected=new, locator=locator,evidence_locator=evidence_locator))
    aliases = c['contig_aliases']
    if len(set(aliases.values())) != len(aliases):
        issue('NONBIJECTIVE_ALIAS',all_ids,aliases,'one-to-one','contract.json#contig_aliases','canonical reference.fa')
    name = lambda v: aliases.get(v,v)
    if {name(k):v for k,v in refs.items()} != erefs:
        issue('REFERENCE_SEQUENCE_CONFLICT',all_ids,'different reference content','official source sequence','reference.fa','SOURCE_LINKS.json; official FASTA SHA256')
    if c['producer_origin'] != c['consumer_origin']:
        issue('DECLARED_COORDINATE_CONFLICT',all_ids,c['consumer_origin'],c['producer_origin'],
              'contract.json#consumer_origin','contract.json#producer_origin; original BAM uses 0-based adapter coordinates')
    if set(c['truth_join_keys']) != set(ec['truth_join_keys']):
        issue('INCOMPLETE_TRUTH_JOIN_KEY',all_ids,c['truth_join_keys'],ec['truth_join_keys'],
              'contract.json#truth_join_keys','sample-specific official control/modified BED mapping')
    if set(c['index_join_keys']) != set(ec['index_join_keys']):
        issue('INCOMPLETE_READ_INDEX_NAMESPACE',all_ids,c['index_join_keys'],ec['index_join_keys'],
              'contract.json#index_join_keys','read_map.tsv; BAM sample/run and UUID')
    def truth_key(row, origin, normalize):
        return (row['sample_id'],normalize(row['contig']),int(row['position'])-origin,
                row['strand'],row['truth_label'],row['source_bed_sha256'],row['source_line'])
    if collections.Counter(truth_key(r,c['truth_origin'],name) for r in truth) != collections.Counter(truth_key(r,0,lambda x:x) for r in etruth):
        issue('OFFICIAL_TRUTH_MAPPING',all_ids,'truth table differs from source mapping','official sample-specific BED records',
              'truth.tsv','contract.json#samples[*].source_files; truth.tsv#source_bed_sha256/source_line in canonical evidence')
    if not read_maps_equal(mapping, emapping):
        issue('UNVERIFIED_READ_MAP',all_ids,'read map differs','verified source-record read map','read_map.tsv','source BAM RG/UUID/ordinal')
    seen = {}
    known_rows = 0
    for row in rows:
        eid = row['evidence_id']
        if eid not in expected:
            issue('UNKNOWN_SOURCE_RECORD',[eid],eid,'verifiable source record',row['_locator'],'no registered source locator')
            continue
        known_rows += 1
        original = expected[eid]
        if eid in seen:
            lhs = {k:v for k,v in row.items() if k not in ('_locator','export_note')}
            rhs = {k:v for k,v in seen[eid].items() if k not in ('_locator','export_note')}
            issue('EXACT_DUPLICATE' if lhs == rhs else 'CONFLICTING_DUPLICATE',[eid],row['_locator'],
                  seen[eid]['_locator'],row['_locator'],f'{evidence.name}/predictions.tsv.gz; {eid}')
        else:
            seen[eid] = row
        location = f"BAM sha256={original['source_bam_sha256']} ordinal={original['bam_record_ordinal']} query_pos0={original['query_pos0']}"
        identity = ('sample_id','run_id','read_id','read_index')
        mismatch = {k:[row[k],original[k]] for k in identity if row[k] != original[k]}
        if mismatch:
            issue('READ_IDENTITY_CONFLICT',[eid],{k:v[0] for k,v in mismatch.items()},
                  {k:v[1] for k,v in mismatch.items()},row['_locator'],location+'; read_map.tsv')
        alignment = {}
        if name(row['contig']) != original['contig']:
            alignment['contig'] = [name(row['contig']),original['contig']]
        if int(row['position'])-c['producer_origin'] != int(original['position']):
            alignment['position'] = [int(row['position'])-c['producer_origin'],int(original['position'])]
        for k in ('strand','query_pos0','read_base','alignment_status','source_bam_sha256','bam_record_ordinal'):
            if row[k] != original[k]:
                alignment[k] = [row[k],original[k]]
        if alignment:
            issue('SOURCE_ALIGNMENT_CONFLICT',[eid],{k:v[0] for k,v in alignment.items()},
                  {k:v[1] for k,v in alignment.items()},row['_locator'],location+'; FASTA/BED exact target')
        if row['modification_code'] not in ('a','28871') or c['modification_code'] not in ('a','28871'):
            issue('MODIFICATION_CODE_CONFLICT',[eid],row['modification_code'],'a or 28871',row['_locator'],location+'; SAMtags modification table')
        if row['score_kind'] != 'read_site_ML_interval':
            issue('SCORE_RESOLUTION_CONFLICT',[eid],row['score_kind'],'read_site_ML_interval',row['_locator'],location+'; endpoint read-site probability')
        if not original['p_mid'] and row['p_mid']:
            issue('MISSING_SCORE_IMPUTED',[eid],row['p_mid'],original['score_status'],row['_locator'],location+'; absent ML score is not exact zero')
        elif any(row[k] != original[k] for k in ('ml_integer','p_low','p_high','p_mid','score_status')):
            issue('SOURCE_SCORE_CONFLICT',[eid],{k:row[k] for k in ('ml_integer','p_mid','score_status')},
                  {k:original[k] for k in ('ml_integer','p_mid','score_status')},row['_locator'],location+'; original MM/ML interval')
    missing = sorted(set(expected)-set(seen))
    if missing or c['expected_opportunities'] != len(expected) or c['population_claim'] != ec['population_claim']:
        issue('POPULATION_SCOPE_CONFLICT',missing,{'observed_unique':len(seen),'declared':c['expected_opportunities'],'claim':c['population_claim']},
              {'expected_unique':len(expected),'claim':ec['population_claim']},'contract.json#population_claim; prediction parts','canonical selected source opportunity ledger, selection before score filtering')
    counts = collections.Counter(x['rule'] for x in issues)
    status = 'UNDETERMINED' if any(k in counts for k in ('UNKNOWN_SOURCE_RECORD','UNVERIFIED_READ_MAP')) else 'REJECT' if issues else 'ACCEPT'
    return dict(method='RMT',status=status,issues=issues,rule_counts=dict(counts),
        input_records=len(rows),unique_source_records=len(seen),expected_opportunities=len(expected),
        explicit_source_probabilities=sum(r['score_status']=='explicit' for r in erows),
        missing_source_probabilities=sum(not r['p_mid'] for r in erows),
        unlocated_input_records=len(rows)-known_rows,unsupported=[],
        scope='source-bound export verification; not discovery of unknown biological mislabelling',
        evidence_file_hashes=links['files'])


def recover(bundle, evidence, output, result):
    rules = set(result['rule_counts'])
    if not rules or not rules <= RECOVERABLE:
        return dict(status='NOT_NEEDED' if not rules else 'REFUSED_NOT_IDENTIFIABLE_OR_OUTSIDE_WHITELIST',
                    rules=sorted(rules),written=False)
    bundle, evidence, output = Path(bundle), Path(evidence), Path(output)
    output.mkdir(parents=True,exist_ok=False)
    c, rows, truth, refs, mapping = load(bundle)
    ec, erows, etruth, erefs, emapping = load(evidence)
    originals = {r['evidence_id']:r for r in erows}
    patches = []
    if 'OFFICIAL_TRUTH_MAPPING' in rules:
        patches.append(dict(file='truth.tsv',old=truth,new=etruth,evidence='official BED SHA256/line and sample mapping'))
        truth = etruth
    if 'DECLARED_COORDINATE_CONFLICT' in rules:
        patches.append(dict(file='contract.json',field='consumer_origin',old=c['consumer_origin'],
                            new=c['producer_origin'],evidence='explicit producer coordinate declaration'))
        c['consumer_origin'] = c['producer_origin']
    unique = {}
    for row in rows:
        eid = row['evidence_id']
        row.pop('_locator')
        original = originals[eid]
        if 'READ_IDENTITY_CONFLICT' in rules:
            for k in ('sample_id','run_id','read_id','read_index'):
                if row[k] != original[k]:
                    patches.append(dict(file='predictions',evidence_id=eid,field=k,old=row[k],new=original[k],
                                        evidence='source BAM SHA256 + ordinal and verified read-index map'))
                    row[k] = original[k]
        if eid in unique:
            if row != unique[eid]:
                raise ValueError('RECOVERY_CONFLICTING_DUPLICATE_REFUSED')
            patches.append(dict(file='predictions',evidence_id=eid,old='one extra identical record',new='one retained',
                                evidence='all fields identical; same native source record and site'))
        else:
            unique[eid] = row
    c['prediction_parts'] = ['predictions.tsv.gz']
    write_table(output/'predictions.tsv.gz',list(unique.values()))
    write_table(output/'truth.tsv',truth)
    write_table(output/'read_map.tsv',mapping)
    write_fasta(output/'reference.fa',refs)
    write_json(output/'contract.json',c)
    post = validate(output,evidence)
    if post['status'] != 'ACCEPT':
        raise ValueError('RECOVERY_POSTVALIDATION_FAILED')
    old_names = load(bundle)[0]['prediction_parts']+['truth.tsv','reference.fa','read_map.tsv','contract.json']
    new_names = c['prediction_parts']+['truth.tsv','reference.fa','read_map.tsv','contract.json']
    detail = dict(status='RECOVERED_SOURCE_DETERMINED',written=True,rules=sorted(rules),patches=patches,
                  before_sha256={n:digest(bundle/n) for n in old_names},after_sha256={n:digest(output/n) for n in new_names},
                  postvalidation_status=post['status'],evidence_sha256=digest(evidence/'SOURCE_LINKS.json'))
    write_json(output/'RECOVERY.json',detail)
    return {k:v for k,v in detail.items() if k != 'patches'}


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for flag in ('bundle','evidence','output'):
        p.add_argument('--'+flag,required=True)
    p.add_argument('--recover')
    a = p.parse_args()
    result = validate(a.bundle,a.evidence)
    if a.recover:
        result['recovery'] = recover(a.bundle,a.evidence,a.recover,result)
    write_json(a.output,result)
