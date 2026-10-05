"""Source adapter and file I/O only; deliberately no candidate/baseline decisions."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, sort_keys=True)
        stream.write('\n')


def table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def write_table(path, rows, fields=None):
    fields = fields or list(rows[0])
    text = io.StringIO(newline='')
    writer = csv.DictWriter(text, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    payload = text.getvalue().encode()
    with Path(path).open('xb') as stream:
        stream.write(gzip.compress(payload, mtime=0) if str(path).endswith('.gz') else payload)


def fasta(path):
    sequences = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            sequences[name] = ''
        elif line.strip():
            sequences[name] += line.strip().upper()
    return sequences


def write_fasta(path, sequences):
    with Path(path).open('x') as stream:
        for name, sequence in sequences.items():
            stream.write(f'>{name}\n{sequence}\n')


def load(path):
    path = Path(path)
    contract = read_json(path/'contract.json')
    rows = []
    for part in contract['prediction_parts']:
        for number, row in enumerate(table(path/part), 2):
            row['_locator'] = f'{part}:{number}'
            rows.append(row)
    return contract, rows, table(path/'truth.tsv'), fasta(path/'reference.fa'), table(path/'read_map.tsv')


def canonical(output, role, canaries):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    observations, truth, mapping, samples = [], [], [], []
    references = None
    for result_dir in canaries:
        result_dir = Path(result_dir).resolve()
        summary = read_json(result_dir/'summary.json')
        if summary['status'] != 'PASS' or summary['role'] != role:
            raise ValueError('CANARY_NOT_PASS_OR_ROLE_CONFLICT')
        source_files = summary['input_files']
        for source in source_files:
            if digest(source['path']) != source['sha256']:
                raise ValueError('CANARY_SOURCE_DRIFT')
        if digest(result_dir/'read_sites.tsv.gz') != summary['record_output_sha256']:
            raise ValueError('CANARY_ROW_DRIFT')
        current_refs = fasta(source_files[1]['path'])
        if references is not None and references != current_refs:
            raise ValueError('PAIR_REFERENCE_DIFFERENCE_REQUIRES_SEPARATE_BUNDLES')
        references = current_refs
        sid = summary['sample_id']
        index_by_read = {}
        for row in table(result_dir/'read_sites.tsv.gz'):
            if row['selected_for_variants'] != '1':
                continue
            read_key = row['run_id'], row['read_id']
            if read_key not in index_by_read:
                index_by_read[read_key] = len(index_by_read)
                mapping.append(dict(sample_id=sid, run_id=row['run_id'],
                                    read_index=index_by_read[read_key], read_id=row['read_id'],
                                    source_bam_sha256=source_files[0]['sha256'],
                                    bam_record_ordinal=row['bam_record_ordinal']))
            item = {k:v for k,v in row.items() if k not in ('truth_label','truth_line','selected_for_variants','ref_pos0')}
            item.update(position=row['ref_pos0'], read_index=str(index_by_read[read_key]),
                        source_bam_sha256=source_files[0]['sha256'],
                        evidence_id='|'.join((sid,row['bam_record_ordinal'],row['contig'],row['ref_pos0'],row['strand'])))
            observations.append(item)
        for line_no, line in enumerate(Path(source_files[2]['path']).read_text().splitlines(), 1):
            if not line or line.startswith('#'):
                continue
            chrom, start, end, code, _, strand = line.split('\t')[:6]
            for pos in range(int(start), int(end)):
                truth.append(dict(sample_id=sid, contig=chrom, position=pos,
                                  strand=strand, truth_label=int(code != '-'),
                                  source_bed_sha256=source_files[2]['sha256'], source_line=line_no))
        samples.append(dict(sample_id=sid, source_files=[dict(path=str(Path(x['path']).relative_to(ROOT)),sha256=x['sha256'],bytes=x['bytes']) for x in source_files],
                            read_groups=summary['read_groups'], original_canary=str(result_dir.relative_to(ROOT)),
                            full_counts=summary['counts'], selected_reads=len(index_by_read)))
    contract = dict(schema_version=1, role=role, bundle_id=output.name,
                    prediction_parts=['predictions.tsv.gz'],
                    producer_origin=0, consumer_origin=0, truth_origin=0,
                    truth_join_keys=['sample_id','contig','position','strand'],
                    index_join_keys=['sample_id','run_id','read_index'],
                    endpoint='READ_SITE_PROBABILITY_EVALUATION', modification_code='a',
                    population_claim='SELECTED_PRIMARY_TRUTH_OVERLAP',
                    selection_rule='first 100 primary mapped unique reads per construct, source BAM order, before score filtering',
                    expected_opportunities=len(observations),
                    contig_aliases={}, samples=samples,
                    independent_training_status='UNKNOWN',
                    redistribution='LOCAL_ONLY_RIGHTS_NOT_CONFIRMED',
                    truth_scope='provider-declared construct condition, not molecular purity/native occupancy')
    write_table(output/'predictions.tsv.gz', observations)
    write_table(output/'truth.tsv', truth)
    write_table(output/'read_map.tsv', mapping)
    write_fasta(output/'reference.fa', references)
    write_json(output/'contract.json', contract)
    write_json(output/'SOURCE_LINKS.json', {'canary_row_sha256':[digest(Path(p)/'read_sites.tsv.gz') for p in canaries],
        'files':{n:digest(output/n) for n in ('contract.json','predictions.tsv.gz','truth.tsv','reference.fa','read_map.tsv')},
        'record_evidence':'Original BAM SHA256 + ordinal + query position; original BED SHA256 + line. Adapter indices are not official m6Anet read indices.'})
    print(json.dumps({'bundle':str(output),'observations':len(observations),'truth_sites':len(truth),'role':role}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',required=True)
    parser.add_argument('--role',required=True)
    parser.add_argument('--canary',action='append',required=True)
    args = parser.parse_args()
    canonical(args.output, args.role, args.canary)
