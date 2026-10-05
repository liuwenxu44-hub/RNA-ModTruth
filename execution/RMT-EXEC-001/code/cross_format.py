"""Inspect actual official m6Anet fixtures, without inference or read-truth claims."""
import argparse
import collections
import csv
import gzip
import json
from pathlib import Path
from bundle import ROOT,digest,write_json


def route(rows,requested):
    columns=set(rows[0])
    if {'read_index','probability_modified','transcript_id','transcript_position'}<=columns:
        level='READ_SITE'
    elif {'n_reads','kmer','mod_ratio','probability_modified','transcript_id','transcript_position'}<=columns:
        level='SITE'
    else:
        raise ValueError('UNSUPPORTED_EXAMPLE_SCHEMA')
    if level!=requested:
        raise ValueError('SCORE_RESOLUTION_CONFLICT:'+level+'->'+requested)
    return level


def run(output):
    files={}
    data={}
    for name in ('data.indiv_proba.csv.gz','data.site_proba.csv.gz'):
        path=ROOT/'.cache/m6anet'/name
        with gzip.open(path,'rt') as stream:
            rows=list(csv.DictReader(stream))
        fields=list(rows[0])
        probabilities=[float(r['probability_modified']) for r in rows]
        if not all(0<=v<=1 for v in probabilities):
            raise ValueError('INVALID_EXAMPLE_PROBABILITY')
        files[name]=dict(sha256=digest(path),bytes=path.stat().st_size,rows=len(rows),columns=fields,
                        probability_range=[min(probabilities),max(probabilities)])
        data[name]=rows
    individual=data['data.indiv_proba.csv.gz']
    site=data['data.site_proba.csv.gz']
    valid_routing=[route(individual,'READ_SITE'),route(site,'SITE')]
    try:
        route(site,'READ_SITE')
        site_misuse='FAILED_TO_REJECT'
    except ValueError as error:
        site_misuse=str(error)
    sid=lambda r:(r['transcript_id'],int(r['transcript_position']))
    readkeys=[sid(r)+(int(r['read_index']),) for r in individual]
    per_site=collections.Counter(sid(r) for r in individual)
    sites={sid(r):r for r in site}
    result=dict(files=files,repository_commit='590ec277cb48d61774f0872395099e466022e810',
        fixture_type='OFFICIAL_SOFTWARE_TEST_OUTPUT_NOT_INDEPENDENT_BIOLOGICAL_VALIDATION',
        duplicate_individual_keys=len(readkeys)-len(set(readkeys)),
        individual_sites=len(per_site),site_rows=len(site),
        individual_sites_without_site_row=len(set(per_site)-set(sites)),
        site_n_reads_matches=sum(int(sites[k]['n_reads'])==n for k,n in per_site.items() if k in sites),
        read_index_uuid_mapping='NOT_AVAILABLE_IN_INSPECTED_OFFICIAL_TEST_TREE',
        read_truth='NOT_AVAILABLE; READ_LEVEL_ACCURACY_NOT_COMPUTED',
        adapter_status='SCHEMA_AND_SITE_LINK_EXECUTED; UUID_AND_TRUTH_JOIN_UNDETERMINED',
        semantics='individual probability_modified and site probability_modified are different resolutions; identical column name is not a conversion rule',
        valid_routing=valid_routing,site_as_read_substitution=site_misuse,
        caller_executed=False)
    write_json(output,result)
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--output',required=True)
    a=p.parse_args()
    run(a.output)
