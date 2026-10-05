"""Bounded self-authored fixtures only; never imports the original repository."""
from contextlib import redirect_stdout
import csv
import io
import json
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from rmt_targeted.bundle import digest, read_json, table
from rmt_targeted.cli import audit, main
from rmt_targeted.example import create_example


class ToolTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(dir=os.environ['RMT_TEST_WORK_ROOT'])
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name) / 'fictional'
        create_example(self.root)
        self.bundle, self.evidence = self.root / 'bundle', self.root / 'evidence'

    def put(self, folder, name, value):
        if name.endswith('.json'):
            (folder / name).write_text(json.dumps(value, allow_nan=False), encoding='utf-8')
        else:
            with (folder / name).open('w', newline='', encoding='utf-8') as stream:
                writer = csv.DictWriter(stream, fieldnames=list(value[0]), delimiter='\t', lineterminator='\n')
                writer.writeheader(); writer.writerows(value)

    def rehash(self):
        links = read_json(self.evidence / 'SOURCE_LINKS.json')
        links['files'] = {name: digest(self.evidence / name) for name in links['files']}
        self.put(self.evidence, 'SOURCE_LINKS.json', links)

    def both(self, status, **kwargs):
        results = []
        for method in ('RMT', 'B2'):
            with self.subTest(method=method):
                result = audit(self.bundle, self.evidence, method, **kwargs)
                self.assertEqual(result['status'], status, result)
                results.append(result)
        self.assertEqual(results[0]['metrics'], results[1]['metrics'])
        return results[0]

    def invoke(self, *args):
        stream = io.StringIO()
        with redirect_stdout(stream):
            code = main(list(args))
        return code, json.loads(stream.getvalue())

    def check_args(self):
        return ['check', '--bundle', str(self.bundle), '--evidence', str(self.evidence)]

    def test_canonical_denominators_and_metrics(self):
        result = self.both('PASS')
        self.assertEqual(result['raw_status'], 'ACCEPT')
        counts = result['metrics']['counts']
        self.assertEqual([counts[k] for k in ('opportunities', 'N', 'O', 'S')], [16, 12, 8, 8])
        self.assertEqual([counts[k] for k in ('excluded_alignment_or_base_scope', 'eligible_missing_probability', 'excluded_by_declared_selection')], [4, 4, 0])
        expected = sum((p - y) ** 2 for p, y in [(16.5 / 256, 0), (240.5 / 256, 0), (32.5 / 256, 1), (224.5 / 256, 1)]) / 4
        self.assertAlmostEqual(result['metrics']['Brier_observed'], expected)

    def test_symmetric_missing_class_not_renormalized(self):
        result = self.both('PASS', selection='symmetric', threshold=.9)
        self.assertEqual(result['metrics']['counts']['S'], 4)
        self.assertIsNone(result['metrics']['Brier_equal_class'])
        self.assertIsNone(result['metrics']['fixed_construct_equal_class_Brier'])
        self.assertEqual(len(result['metrics']['missing_required_construct_classes']), 2)

    def test_one_sided_is_only_conditional(self):
        result = self.both('PASS', selection='one-sided', threshold=.9)
        self.assertEqual(result['metrics']['counts']['S'], 2)
        self.assertFalse(result['metrics']['full_population_loss_identified'])

    def test_empty_selection_is_insufficient(self):
        result = self.both('INSUFFICIENT', selection='symmetric', threshold=1.)
        self.assertEqual(result['audit_status'], 'PASS')
        self.assertIsNone(result['metrics']['Brier_observed'])

    def test_row_and_map_order_invariance(self):
        for name in ('predictions.tsv', 'read_map.tsv', 'truth.tsv'):
            self.put(self.bundle, name, list(reversed(table(self.bundle / name))))
        self.both('PASS')

    def test_exact_candidate_duplicate_rejected(self):
        rows = table(self.bundle / 'predictions.tsv'); rows.append(dict(rows[0]))
        self.put(self.bundle, 'predictions.tsv', rows)
        self.both('REJECT')

    def test_conflicting_candidate_duplicate_rejected(self):
        rows = table(self.bundle / 'predictions.tsv'); extra = dict(rows[0]); extra['read_id'] = 'fictional_conflict'; rows.append(extra)
        self.put(self.bundle, 'predictions.tsv', rows)
        self.both('REJECT')

    def test_candidate_map_duplicate_is_insufficient(self):
        rows = table(self.bundle / 'read_map.tsv'); rows.append(dict(rows[0]))
        self.put(self.bundle, 'read_map.tsv', rows)
        self.both('INSUFFICIENT')

    def test_candidate_map_conflict_is_insufficient(self):
        rows = table(self.bundle / 'read_map.tsv'); rows[0]['read_id'] = 'fictional_conflict'
        self.put(self.bundle, 'read_map.tsv', rows)
        self.both('INSUFFICIENT')

    def test_source_evidence_duplicate_never_deduplicated(self):
        rows = table(self.evidence / 'predictions.tsv'); rows.append(dict(rows[0]))
        self.put(self.evidence, 'predictions.tsv', rows); self.rehash()
        result = self.both('INSUFFICIENT')
        self.assertEqual(result['issues'][0]['rule'], 'DUPLICATE_SOURCE_EVIDENCE_ID')

    def test_source_evidence_conflict_never_deduplicated(self):
        rows = table(self.evidence / 'predictions.tsv'); extra = dict(rows[0]); extra['read_id'] = 'fictional_conflict'; rows.append(extra)
        self.put(self.evidence, 'predictions.tsv', rows); self.rehash()
        result = self.both('INSUFFICIENT')
        self.assertEqual(result['issues'][0]['rule'], 'CONFLICTING_SOURCE_EVIDENCE_ID')

    def test_source_composite_collision_with_different_id(self):
        rows = table(self.evidence / 'predictions.tsv'); extra = dict(rows[0]); extra['evidence_id'] = 'fictional_other_id'; rows.append(extra)
        self.put(self.evidence, 'predictions.tsv', rows); self.rehash()
        result = self.both('INSUFFICIENT')
        self.assertEqual(result['issues'][0]['rule'], 'CONFLICTING_SOURCE_COMPOSITE_KEY')

    def test_conflicting_source_truth(self):
        rows = table(self.evidence / 'truth.tsv'); extra = dict(rows[0]); extra['truth_label'] = '1'; rows.append(extra)
        self.put(self.evidence, 'truth.tsv', rows); self.rehash()
        self.both('INSUFFICIENT')

    def test_conflicting_candidate_truth(self):
        rows = table(self.bundle / 'truth.tsv'); extra = dict(rows[0]); extra['truth_label'] = '1'; rows.append(extra)
        self.put(self.bundle, 'truth.tsv', rows)
        self.both('REJECT')

    def test_source_map_duplicate(self):
        rows = table(self.evidence / 'read_map.tsv'); rows.append(dict(rows[0]))
        self.put(self.evidence, 'read_map.tsv', rows); self.rehash()
        self.both('INSUFFICIENT')

    def test_unclosed_source_hash_inventory(self):
        links = read_json(self.evidence / 'SOURCE_LINKS.json'); links['files'].pop('truth.tsv')
        self.put(self.evidence, 'SOURCE_LINKS.json', links)
        self.both('INSUFFICIENT')

    def test_source_hash_drift(self):
        with (self.evidence / 'truth.tsv').open('a') as stream: stream.write('\n')
        self.both('INSUFFICIENT')

    def test_duplicate_json_key(self):
        path = self.evidence / 'SOURCE_LINKS.json'
        path.write_text('{"files":{},"files":{}}', encoding='utf-8')
        self.both('INSUFFICIENT')

    def test_unknown_candidate_source_is_insufficient(self):
        rows = table(self.bundle / 'predictions.tsv'); rows[0]['evidence_id'] = 'fictional_unknown'
        self.put(self.bundle, 'predictions.tsv', rows)
        self.both('INSUFFICIENT')

    def test_truncation_does_not_become_legal_by_metric_selection(self):
        rows = table(self.bundle / 'predictions.tsv')
        self.put(self.bundle, 'predictions.tsv', rows[1:])
        self.both('REJECT', selection='symmetric', threshold=.9)

    def test_malformed_candidate_schema(self):
        rows = table(self.bundle / 'predictions.tsv')
        for row in rows: del row['run_id']
        self.put(self.bundle, 'predictions.tsv', rows)
        self.both('REJECT')

    def test_no_scores_has_null_metrics(self):
        for folder in (self.bundle, self.evidence):
            rows = table(folder / 'predictions.tsv')
            for row in rows:
                for field in ('p_mid', 'p_low', 'p_high', 'ml_integer'): row[field] = ''
                if row['score_status'] == 'explicit': row['score_status'] = 'fictional_missing_unspecified'
            self.put(folder, 'predictions.tsv', rows)
        self.rehash()
        result = self.both('INSUFFICIENT')
        self.assertEqual(result['audit_status'], 'PASS')
        self.assertEqual(result['metrics']['counts']['O'], 0)
        self.assertIsNone(result['metrics']['Brier_observed'])

    def test_default_and_strict_rejection_exit(self):
        rows = table(self.bundle / 'predictions.tsv'); rows.append(dict(rows[0]))
        self.put(self.bundle, 'predictions.tsv', rows)
        code, data = self.invoke(*self.check_args())
        self.assertEqual((code, data['status']), (0, 'REJECT'))
        code, data = self.invoke(*self.check_args(), '--strict-exit')
        self.assertEqual((code, data['status']), (1, 'REJECT'))

    def test_pass_and_insufficient_exit(self):
        code, data = self.invoke(*self.check_args(), '--strict-exit')
        self.assertEqual((code, data['status']), (0, 'PASS'))
        (self.evidence / 'SOURCE_LINKS.json').unlink()
        code, data = self.invoke(*self.check_args())
        self.assertEqual((code, data['status']), (0, 'INSUFFICIENT'))
        code, data = self.invoke(*self.check_args(), '--strict-exit')
        self.assertEqual((code, data['status']), (2, 'INSUFFICIENT'))

    def test_usage_exit(self):
        for extra in (['--selection', 'symmetric'], ['--selection', 'symmetric', '--threshold', 'nan'], ['--bad-flag']):
            code, data = self.invoke(*self.check_args(), *extra)
            self.assertEqual((code, data['status']), (64, 'USAGE_ERROR'))

    def test_output_is_never_overwritten(self):
        target = self.root / 'existing.json'; target.write_text('user-owned marker', encoding='utf-8')
        code, data = self.invoke(*self.check_args(), '--output', str(target))
        self.assertEqual((code, data['status']), (3, 'EXECUTION_FAILURE'))
        self.assertEqual(target.read_text(), 'user-owned marker')

    def test_runtime_valueerror_is_not_usage_error(self):
        with patch('rmt_targeted.cli.audit', side_effect=ValueError('FICTIONAL_INJECTED_RUNTIME_FAULT')):
            code, data = self.invoke(*self.check_args())
        self.assertEqual((code, data['status']), (3, 'EXECUTION_FAILURE'))

    def test_unsafe_input_path(self):
        contract = read_json(self.bundle / 'contract.json'); contract['prediction_parts'] = ['../escape.tsv']
        self.put(self.bundle, 'contract.json', contract)
        self.both('REJECT')

    def test_no_truth_join_source(self):
        rows = table(self.evidence / 'truth.tsv'); self.put(self.evidence, 'truth.tsv', rows[1:]); self.rehash()
        self.both('INSUFFICIENT')


if __name__ == '__main__':
    unittest.main(verbosity=2)
