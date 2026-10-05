"""Shared evidence/schema support; no method receives extra source information."""
import collections
import math
from pathlib import Path
import re

from .bundle import digest, load, read_json

FIELDS = ('evidence_id', 'sample_id', 'run_id', 'read_id', 'read_index', 'contig',
          'position', 'strand', 'query_pos0', 'read_base', 'alignment_status',
          'source_bam_sha256', 'bam_record_ordinal', 'modification_code', 'score_kind',
          'ml_integer', 'p_low', 'p_high', 'p_mid', 'score_status')
MAP_FIELDS = ('sample_id', 'run_id', 'read_index', 'read_id', 'source_bam_sha256', 'bam_record_ordinal')
TRUTH_FIELDS = ('sample_id', 'contig', 'position', 'strand', 'truth_label', 'source_bed_sha256', 'source_line')
SOURCE_KEY = ('sample_id', 'run_id', 'source_bam_sha256', 'bam_record_ordinal',
              'query_pos0', 'contig', 'position', 'strand')
READ_KEY = ('sample_id', 'run_id', 'read_index')
TRUTH_KEY = ('sample_id', 'contig', 'position', 'strand')
REQUIRED_CONTRACT = ('schema_version', 'prediction_parts', 'producer_origin', 'consumer_origin',
                     'truth_origin', 'truth_join_keys', 'index_join_keys', 'endpoint',
                     'modification_code', 'population_claim', 'expected_opportunities', 'contig_aliases')


class InputProblem(Exception):
    def __init__(self, code, source, detail):
        self.code, self.source, self.detail = code, source, detail
        super().__init__(code + ': ' + detail)


def key(row, fields):
    """Tuple identity, never delimiter concatenation or row-position identity."""
    return tuple(row[field] for field in fields)


def _unique(rows, fields, kind, fail):
    indexed = {}
    for row in rows:
        identity = key(row, fields)
        if identity in indexed:
            same = {k: v for k, v in row.items() if k != '_locator'} == indexed[identity]
            fail(('DUPLICATE_' if same else 'CONFLICTING_') + kind, repr(identity))
        indexed[identity] = {k: v for k, v in row.items() if k != '_locator'}


def inspect(path, source=False):
    path = Path(path).resolve()
    def fail(code, detail):
        raise InputProblem(code, source, str(detail))
    try:
        c = read_json(path / 'contract.json')
        if not isinstance(c, dict) or any(k not in c for k in REQUIRED_CONTRACT):
            fail('INCOMPLETE_CONTRACT', 'Required version-1 canonical fields missing')
        if c['schema_version'] != 1 or c['endpoint'] != 'READ_SITE_PROBABILITY_EVALUATION':
            fail('UNSUPPORTED_CONTRACT', 'Only the existing schema-v1 read-site endpoint is supported')
        parts = c['prediction_parts']
        if not isinstance(parts, list) or not parts or len(parts) != len(set(parts)):
            fail('INVALID_PREDICTION_PARTS', 'Nonempty unique relative filenames required')
        names = parts + ['truth.tsv', 'reference.fa', 'read_map.tsv', 'contract.json']
        if len(names) != len(set(names)):
            fail('OVERLAPPING_INPUT_ROLES', 'Prediction parts may not overwrite other bundle roles')
        for name in names:
            if not isinstance(name, str) or Path(name).is_absolute() or '..' in Path(name).parts:
                fail('UNSAFE_BUNDLE_PATH', name)
            if not (path / name).resolve().is_relative_to(path):
                fail('BUNDLE_PATH_ESCAPE', name)
        for field in ('producer_origin', 'consumer_origin', 'truth_origin'):
            if type(c[field]) is not int or c[field] not in (0, 1):
                fail('INVALID_COORDINATE_ORIGIN', field)
        if source and any(c[k] != 0 for k in ('producer_origin', 'consumer_origin', 'truth_origin')):
            fail('UNSUPPORTED_SOURCE_COORDINATES', 'Existing source ledgers are zero-based')
        if type(c['expected_opportunities']) is not int or c['expected_opportunities'] < 0:
            fail('INVALID_DECLARED_DENOMINATOR', 'Nonnegative integer required')
        if not isinstance(c['population_claim'], str) or not c['population_claim']:
            fail('MISSING_POPULATION_SCOPE', 'An explicit population claim is required')
        if not isinstance(c['contig_aliases'], dict) or any(not isinstance(x, str) for pair in c['contig_aliases'].items() for x in pair):
            fail('INVALID_CONTIG_ALIASES', 'String-to-string dictionary required')
        for field, allowed in [('truth_join_keys', set(TRUTH_FIELDS)), ('index_join_keys', set(MAP_FIELDS))]:
            vals = c[field]
            if not isinstance(vals, list) or not vals or len(vals) != len(set(vals)) or not set(vals) <= allowed:
                fail('INVALID_JOIN_KEY_DECLARATION', field)
        if source:
            links = read_json(path / 'SOURCE_LINKS.json')
            if not isinstance(links, dict) or not isinstance(links.get('files'), dict):
                fail('INCOMPLETE_SOURCE_HASH_MANIFEST', 'SOURCE_LINKS.json files object required')
            if set(links['files']) != set(names):
                fail('SOURCE_HASH_INVENTORY_NOT_CLOSED', 'Manifest must hash all and only declared bundle files')
            for name, sha in links['files'].items():
                if not isinstance(sha, str) or not re.fullmatch('[0-9a-f]{64}', sha) or digest(path / name) != sha:
                    fail('SOURCE_LEDGER_HASH_MISMATCH', name)
        data = load(path)
        _, rows, truth, refs, mapping = data
        for records, fields, title in [(rows, FIELDS, 'predictions'), (truth, TRUTH_FIELDS, 'truth'), (mapping, MAP_FIELDS, 'read_map')]:
            for row in records:
                if None in row or any(field not in row or row[field] is None for field in fields):
                    fail('MALFORMED_TABLE_SCHEMA', title)
                nullable = {'query_pos0', 'read_base', 'ml_integer', 'p_low', 'p_high', 'p_mid'}
                if any(row[field] == '' for field in fields if field not in nullable):
                    fail('MISSING_REQUIRED_IDENTITY', title)
        if source and (not rows or not truth or not mapping or not refs):
            fail('EMPTY_SOURCE_LEDGER', 'A nonempty complete source ledger is required')
        _unique(truth, TRUTH_KEY, 'TRUTH_KEY', fail)
        # Known candidate map defects remain UNDETERMINED in both immutable checker copies.
        if source:
            _unique(rows, ('evidence_id',), 'SOURCE_EVIDENCE_ID', fail)
            _unique(rows, SOURCE_KEY, 'SOURCE_COMPOSITE_KEY', fail)
            _unique(mapping, READ_KEY, 'SOURCE_READ_MAP_KEY', fail)
        for row in truth:
            if row['truth_label'] not in ('0', '1') or row['strand'] not in ('+', '-'):
                fail('INVALID_TRUTH_VALUE', repr(key(row, TRUTH_KEY)))
            int(row['position']); int(row['source_line'])
        for row in rows:
            int(row['position']); int(row['bam_record_ordinal']); int(row['read_index'])
            if row['query_pos0']:
                int(row['query_pos0'])
            if row['strand'] not in ('+', '-'):
                fail('INVALID_STRAND', row['evidence_id'])
            if row['p_mid']:
                for field in ('p_low', 'p_high', 'p_mid'):
                    v = float(row[field])
                    if not math.isfinite(v) or not 0 <= v <= 1:
                        fail('INVALID_PROBABILITY', row['evidence_id'])
            if row['ml_integer'] and not 0 <= int(row['ml_integer']) <= 255:
                fail('INVALID_ML_INTEGER', row['evidence_id'])
        if source:
            if c['truth_join_keys'] != list(TRUTH_KEY) or c['index_join_keys'] != list(READ_KEY):
                fail('UNSUPPORTED_SOURCE_KEYS', 'Full composite source keys required')
            if c['expected_opportunities'] != len(rows):
                fail('SOURCE_DENOMINATOR_CONFLICT', 'Declared opportunities do not equal unique source rows')
            labels = {key(r, TRUTH_KEY): r for r in truth}
            maps = {key(r, READ_KEY): r for r in mapping}
            for row in rows:
                if key(row, TRUTH_KEY) not in labels:
                    fail('SOURCE_TRUTH_JOIN_NOT_UNIQUE', row['evidence_id'])
                mapped = maps.get(key(row, READ_KEY))
                if mapped is None or any(mapped[k] != row[k] for k in MAP_FIELDS):
                    fail('SOURCE_READ_MAP_CONFLICT', row['evidence_id'])
                if row['contig'] not in refs or not 0 <= int(row['position']) < len(refs[row['contig']]):
                    fail('SOURCE_REFERENCE_COORDINATE_INVALID', row['evidence_id'])
                if row['p_mid']:
                    m = int(row['ml_integer'])
                    if (row['score_status'] != 'explicit' or row['alignment_status'] != 'MATCHED_READ_SITE'
                            or row['score_kind'] != 'read_site_ML_interval' or row['modification_code'] not in ('a', '28871')
                            or float(row['p_mid']) != (m + .5) / 256 or float(row['p_low']) != m / 256
                            or float(row['p_high']) != (m + 1) / 256):
                        fail('SOURCE_SCORE_SEMANTICS_CONFLICT', row['evidence_id'])
                elif row['score_status'] == 'explicit' or any(row[k] for k in ('ml_integer', 'p_low', 'p_high')):
                    fail('SOURCE_MISSING_STATE_CONFLICT', row['evidence_id'])
        return data
    except InputProblem:
        raise
    except (FileNotFoundError, ValueError, KeyError, TypeError, IndexError) as exc:
        fail('SOURCE_INPUT_INSUFFICIENT' if source else 'MALFORMED_BUNDLE', type(exc).__name__ + ': ' + str(exc))
