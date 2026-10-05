"""Entirely self-authored fictional fixture; contains no original study records."""
import hashlib
import json
from pathlib import Path

from .bundle import digest, write_fasta, write_json, write_table


def create_example(output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    refs = {'fictional_construct_A': 'CAAAAAAC', 'fictional_construct_B': 'GAAAAAAG'}
    observations, truth, mapping = [], [], []
    for y in (0, 1):
        sample = 'fictional_control' if y == 0 else 'fictional_modified'
        source = hashlib.sha256(('FICTIONAL_NOT_A_REAL_BAM:' + sample).encode()).hexdigest()
        bed = hashlib.sha256(('FICTIONAL_NOT_A_REAL_BED:' + sample).encode()).hexdigest()
        ordinal = 0
        for construct in refs:
            truth.append(dict(sample_id=sample, contig=construct, position='2', strand='+', truth_label=str(y),
                              source_bed_sha256=bed, source_line=str(len(truth) + 1)))
            for j in range(4):
                ordinal += 1
                identity = dict(sample_id=sample, run_id='fictional_run', read_index=str(ordinal),
                                read_id=f'fictional_read_{y}_{ordinal}', source_bam_sha256=source,
                                bam_record_ordinal=str(ordinal))
                mapping.append(dict(identity))
                m = ([16, 240] if y == 0 else [32, 224])[j] if j < 2 else None
                native_key = [sample, identity['run_id'], source, str(ordinal), '2', construct, '2', '+']
                eid = 'fictional:' + hashlib.sha256(json.dumps(native_key, separators=(',', ':')).encode()).hexdigest()
                observations.append(dict(identity, evidence_id=eid, contig=construct, position='2', strand='+',
                    query_pos0='2', read_base='A' if j < 3 else 'C',
                    alignment_status='MATCHED_READ_SITE' if j < 3 else 'CANONICAL_BASE_MISMATCH',
                    modification_code='a', score_kind='read_site_ML_interval',
                    ml_integer=str(m) if m is not None else '', p_low=str(m / 256) if m is not None else '',
                    p_high=str((m + 1) / 256) if m is not None else '',
                    p_mid=str((m + .5) / 256) if m is not None else '',
                    score_status='explicit' if m is not None else 'fictional_missing_unspecified' if j == 2 else 'fictional_base_mismatch'))
    contract = dict(schema_version=1, bundle_id='self_authored_fictional_example', role='FICTIONAL_TEST_ONLY',
        prediction_parts=['predictions.tsv'], producer_origin=0, consumer_origin=0, truth_origin=0,
        truth_join_keys=['sample_id', 'contig', 'position', 'strand'],
        index_join_keys=['sample_id', 'run_id', 'read_index'], endpoint='READ_SITE_PROBABILITY_EVALUATION',
        modification_code='a', population_claim='FICTIONAL_SAVED_TRUTH_OVERLAP', expected_opportunities=len(observations),
        contig_aliases={}, truth_scope='Invented labels; no biological validation', samples=[],
        independent_training_status='NOT_APPLICABLE_FICTIONAL', redistribution='SELF_AUTHORED_FIXTURE_NOT_THIRD_PARTY_DATA',
        source_locator_notice='BAM/BED SHA-shaped tokens are fictional identifiers, not hashes of existing native files')
    for directory in ('bundle', 'evidence'):
        folder = output / directory
        folder.mkdir()
        write_json(folder / 'contract.json', contract)
        write_table(folder / 'predictions.tsv', observations)
        write_table(folder / 'truth.tsv', truth)
        write_table(folder / 'read_map.tsv', mapping)
        write_fasta(folder / 'reference.fa', refs)
        if directory == 'evidence':
            write_json(folder / 'SOURCE_LINKS.json', dict(
                files={name: digest(folder / name) for name in ('contract.json', 'predictions.tsv', 'truth.tsv', 'read_map.tsv', 'reference.fa')},
                record_evidence='Entirely invented fixture. These hashes certify local byte consistency only, not native provenance.'))
    write_json(output / 'FICTIONAL_EXAMPLE.json', dict(
        status='FICTIONAL_TEST_ONLY', original_study_records=0, original_UUIDs=0,
        native_BAM_BED_files=0, description='Two invented sequences, invented read labels and deliberately chosen integer scores.',
        license='Project software: MIT (repository-root LICENSE); this example contains no third-party study records.'))
    return dict(status='PASS', operation='CREATE_FICTIONAL_EXAMPLE', output=str(output.resolve()), fictional=True)
