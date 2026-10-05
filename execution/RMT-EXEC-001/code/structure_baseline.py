"""B0: explicit structural/schema-equivalent validation, not endpoint inference."""
import argparse
import math
from pathlib import Path
from bundle import load, write_json


def check(path):
    issues=[]
    try:
        c,rows,truth,refs,mapping=load(path)
        required={'schema_version':int,'role':str,'prediction_parts':list,
                  'producer_origin':int,'consumer_origin':int,'truth_origin':int,
                  'truth_join_keys':list,'index_join_keys':list,'samples':list,
                  'endpoint':str,'modification_code':str,'population_claim':str,
                  'expected_opportunities':int,'contig_aliases':dict}
        for name,kind in required.items():
            if name not in c or type(c[name]) is not kind:
                raise ValueError('CONTRACT_TYPE:'+name)
        if c['schema_version'] != 1 or any(c[x] not in (0,1) for x in ('producer_origin','consumer_origin','truth_origin')):
            raise ValueError('SCHEMA_VERSION_OR_COORDINATE_ENUM')
        if not c['samples'] or not refs or not rows or not truth:
            raise ValueError('EMPTY_REQUIRED_TABLE')
        need=('evidence_id','sample_id','run_id','read_id','read_index','contig','position',
              'strand','query_pos0','score_kind','score_status','modification_code',
              'ml_integer','p_low','p_high','p_mid','source_bam_sha256','bam_record_ordinal')
        for row in rows:
            if any(k not in row for k in need):
                raise ValueError('PREDICTION_REQUIRED_COLUMN')
            if row['strand'] not in ('+','-'):
                raise ValueError('STRAND_ENUM')
            for k in ('position','read_index','bam_record_ordinal'):
                if int(row[k])<0:
                    raise ValueError('NEGATIVE_INTEGER:'+k)
            if row['ml_integer'] and not 0<=int(row['ml_integer'])<=255:
                raise ValueError('ML_BYTE_RANGE')
            for k in ('p_low','p_high','p_mid'):
                if row[k] and (not math.isfinite(float(row[k])) or not 0<=float(row[k])<=1):
                    raise ValueError('PROBABILITY_RANGE:'+k)
        for row in truth:
            if row['truth_label'] not in ('0','1') or row['strand'] not in ('+','-') or int(row['position'])<0:
                raise ValueError('TRUTH_STRUCTURAL_DOMAIN')
    except Exception as exc:
        issues.append(dict(rule='STRUCTURAL_INVALID',evidence_ids=[],locator=str(path),evidence_locator=str(exc)))
    return dict(method='B0',status='REJECT' if issues else 'ACCEPT',issues=issues,
        rule_counts={'STRUCTURAL_INVALID':len(issues)} if issues else {},
        input_records=len(rows) if 'rows' in locals() else None,
        unsupported=['scientific endpoint','truth/probability semantics','cross-file scientific identity','recovery'],
        scope='STRUCTURAL_ONLY_ENDPOINT_UNASSESSED; ACCEPT means schema/readability only')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('bundle','evidence','output'):
        p.add_argument('--'+k,required=True)
    a=p.parse_args()
    write_json(a.output,check(a.bundle))
