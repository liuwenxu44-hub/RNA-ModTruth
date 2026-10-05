# Versioned validation-only copy; original source remains immutable.
# Parent: methods/RMT-METHODS-001/code/checkers_v2/simple_baseline.py
"""B2: independent SQL endpoint checker, with independent allowed recovery.

Imports only shared serialization. No RMT decision/recovery functions are imported.
Uses the same source ledger, facts, budget and development cases as RMT.
"""
import argparse
import collections
import sqlite3
import json
from pathlib import Path
from .bundle import load, read_json, digest, write_json, write_table, write_fasta

FIELDS = ['evidence_id','sample_id','run_id','read_id','read_index','contig','position',
          'strand','query_pos0','read_base','alignment_status','source_bam_sha256',
          'bam_record_ordinal','modification_code','score_kind','ml_integer','p_low','p_high',
          'p_mid','score_status']


def read_maps_equal(observed, source):
    """Independent SQL implementation of orderless source-map comparison."""
    sql = sqlite3.connect(':memory:')
    try:
        for name, rows in [('actual_map', observed), ('source_map', source)]:
            sql.execute('CREATE TABLE '+name+' (sample TEXT, run TEXT, idx TEXT, payload TEXT)')
            for row in rows:
                if any(not row.get(k, '') for k in ('sample_id', 'run_id', 'read_index')):
                    return False
                sql.execute('INSERT INTO '+name+' VALUES (?,?,?,?)',
                    (row['sample_id'], row['run_id'], row['read_index'], json.dumps(row,sort_keys=True,separators=(',',':'))))
            if sql.execute('SELECT 1 FROM '+name+' GROUP BY sample,run,idx HAVING COUNT(*)>1 LIMIT 1').fetchone():
                return False
        for a,b in [('actual_map','source_map'),('source_map','actual_map')]:
            if sql.execute('SELECT * FROM '+a+' EXCEPT SELECT * FROM '+b).fetchone():
                return False
        return True
    finally:
        sql.close()


def check(bundle, evidence):
    c, data, truth, refs, readmap = load(bundle)
    e, parent, goldtruth, goldrefs, goldmap = load(evidence)
    source_links = read_json(Path(evidence)/'SOURCE_LINKS.json')
    if not all(digest(Path(evidence)/k) == v for k,v in source_links['files'].items()):
        raise ValueError('SOURCE_LEDGER_HASH_MISMATCH')
    db = sqlite3.connect(':memory:')
    for tab in ('observed','expected'):
        db.execute('CREATE TABLE '+tab+' ('+','.join(k+' TEXT NOT NULL' for k in FIELDS)+', locator TEXT)')
    def values(records, contract):
        aliases = contract['contig_aliases']
        for row in records:
            z = dict(row)
            z['contig'] = aliases.get(row['contig'],row['contig'])
            z['position'] = str(int(row['position'])-contract['producer_origin'])
            yield tuple(z[k] for k in FIELDS)+(row['_locator'],)
    marks = ','.join('?' for _ in range(len(FIELDS)+1))
    db.executemany('INSERT INTO observed VALUES ('+marks+')',values(data,c))
    db.executemany('INSERT INTO expected VALUES ('+marks+')',values(parent,e))
    db.execute('CREATE INDEX expected_id ON expected(evidence_id)')
    db.execute('CREATE INDEX observed_id ON observed(evidence_id)')
    issues = []
    allids = [x['evidence_id'] for x in parent]
    def add(code, ids, location, reason):
        ids = sorted(set(ids))
        if ids:
            issues.append(dict(rule=code,evidence_ids=ids,locator=location,
                               evidence_locator=reason))
    def query(code, condition, reason):
        found = db.execute('SELECT DISTINCT o.evidence_id FROM observed o JOIN expected e USING(evidence_id) WHERE '+condition).fetchall()
        add(code,[r[0] for r in found],'prediction_parts; source evidence IDs',reason)
    norm = lambda x: c['contig_aliases'].get(x,x)
    if len(set(c['contig_aliases'].values())) != len(c['contig_aliases']):
        add('ALIAS_CONFLICT',allids,'contract.json#contig_aliases','documented bijections only')
    if {norm(k):v for k,v in refs.items()} != goldrefs:
        add('FASTA_CONFLICT',allids,'reference.fa','official reference sequence SHA256')
    if c['consumer_origin'] != c['producer_origin']:
        add('COORDINATE_DECLARATION',allids,'contract.json#consumer_origin','explicit producer-origin declaration')
    if sorted(c['truth_join_keys']) != sorted(e['truth_join_keys']):
        add('JOIN_KEYS',allids,'contract.json#truth_join_keys','sample-specific official BED mapping requires sample/contig/position/strand')
    if sorted(c['index_join_keys']) != sorted(e['index_join_keys']):
        add('INDEX_NAMESPACE',allids,'contract.json#index_join_keys','read index requires sample and run namespace')
    def tkey(r, origin, name):
        return (r['sample_id'],name(r['contig']),int(r['position'])-origin,r['strand'],r['truth_label'],r['source_bed_sha256'],r['source_line'])
    if collections.Counter(tkey(r,c['truth_origin'],norm) for r in truth) != collections.Counter(tkey(r,0,lambda x:x) for r in goldtruth):
        add('SAMPLE_BED_MAPPING',allids,'truth.tsv','source BED SHA256/line and official sample association')
    map_equal = read_maps_equal(readmap, goldmap)
    if not map_equal:
        add('READ_MAP_UNVERIFIED',allids,'read_map.tsv','source BAM header/read UUID and original ordinal')
    unknown = db.execute('SELECT DISTINCT evidence_id FROM observed EXCEPT SELECT evidence_id FROM expected').fetchall()
    add('UNLOCATED_RECORD',[x[0] for x in unknown],'predictions','no source attestation for new records')
    query('IDENTITY', ' OR '.join('o.'+k+'<>e.'+k for k in ('sample_id','run_id','read_id','read_index')),
          'native BAM source-record locator and read_map.tsv')
    query('ALIGNMENT', ' OR '.join('o.'+k+'<>e.'+k for k in ('contig','position','strand','query_pos0','read_base','alignment_status','source_bam_sha256','bam_record_ordinal')),
          'native BAM CIGAR/query position; official FASTA and BED')
    query('MOD_CODE', "o.modification_code NOT IN ('a','28871')",'SAMtags code a/28871 is m6A, not 17596')
    if c['modification_code'] not in ('a','28871'):
        add('MOD_CODE',allids,'contract.json#modification_code','requested m6A endpoint')
    query('SCORE_KIND', "o.score_kind<>'read_site_ML_interval'",'read-site endpoint cannot consume site-aggregate probabilities')
    query('MISSING_IMPUTED', "e.p_mid='' AND o.p_mid<>''",'absence of native ML probability is not zero')
    query('SCORE_VALUE', "NOT (e.p_mid='' AND o.p_mid<>'') AND ("+' OR '.join('o.'+k+'<>e.'+k for k in ('ml_integer','p_low','p_high','p_mid','score_status'))+')',
          'original MM/ML integer, interval and missing-state source record')
    counts = collections.Counter(r['evidence_id'] for r in data)
    for eid,n in counts.items():
        if n > 1:
            repeated = [tuple(r.get(k,'') for k in FIELDS) for r in data if r['evidence_id']==eid]
            add('EXACT_REPEAT' if len(set(repeated))==1 else 'CONFLICTING_REPEAT',[eid],
                'prediction_parts','same original record and site; all scientific fields compared')
    missing = [x[0] for x in db.execute('SELECT evidence_id FROM expected EXCEPT SELECT evidence_id FROM observed')]
    if missing or c['expected_opportunities'] != len(parent) or c['population_claim'] != e['population_claim']:
        add('DENOMINATOR',missing or allids,'contract.json#population_claim; prediction_parts',
            'source opportunities selected before filtering, including absent score states')
    db.close()
    status = 'UNDETERMINED' if unknown or not map_equal else 'REJECT' if issues else 'ACCEPT'
    return dict(method='B2',status=status,issues=issues,unsupported=[],
        input_records=len(data),unique_source_records=len(counts),expected_opportunities=len(parent),
        explicit_source_probabilities=sum(r['score_status']=='explicit' for r in parent),
        missing_source_probabilities=sum(not r['p_mid'] for r in parent),
        rule_counts=dict(collections.Counter(x['rule'] for x in issues)),
        scope='independent SQL source-bound endpoint checker; same facts and recovery whitelist as RMT')


