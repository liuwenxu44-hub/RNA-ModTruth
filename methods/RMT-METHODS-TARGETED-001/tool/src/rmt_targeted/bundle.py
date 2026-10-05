"""Path-independent serialization shared by both validation implementations."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def read_json(path):
    def unique_pairs(pairs):
        out = {}
        for key, value in pairs:
            if key in out:
                raise ValueError('DUPLICATE_JSON_KEY: ' + key)
            out[key] = value
        return out
    return json.loads(Path(path).read_text(encoding='utf-8'), object_pairs_hook=unique_pairs)


def write_json(path, value):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, sort_keys=True, allow_nan=False)
        stream.write('\n')


def table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', encoding='utf-8', newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        fields = reader.fieldnames or []
        if not fields or len(fields) != len(set(fields)):
            raise ValueError('MISSING_OR_DUPLICATE_TABLE_HEADER')
        return list(reader)


def write_table(path, rows, fields=None):
    fields = fields or list(rows[0])
    text = io.StringIO(newline='')
    writer = csv.DictWriter(text, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    payload = text.getvalue().encode('utf-8')
    with Path(path).open('xb') as stream:
        stream.write(gzip.compress(payload, mtime=0) if str(path).endswith('.gz') else payload)


def fasta(path):
    sequences = {}
    name = None
    for line in Path(path).read_text(encoding='utf-8').splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in sequences:
                raise ValueError('DUPLICATE_FASTA_IDENTIFIER')
            sequences[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('FASTA_SEQUENCE_WITHOUT_HEADER')
            sequences[name] += line.strip().upper()
    return sequences


def write_fasta(path, sequences):
    with Path(path).open('x', encoding='utf-8') as stream:
        for name, sequence in sequences.items():
            stream.write(f'>{name}\n{sequence}\n')


def load(path):
    path = Path(path)
    contract = read_json(path / 'contract.json')
    rows = []
    for part in contract['prediction_parts']:
        for number, row in enumerate(table(path / part), 2):
            row['_locator'] = f'{part}:{number}'
            rows.append(row)
    return contract, rows, table(path / 'truth.tsv'), fasta(path / 'reference.fa'), table(path / 'read_map.tsv')
