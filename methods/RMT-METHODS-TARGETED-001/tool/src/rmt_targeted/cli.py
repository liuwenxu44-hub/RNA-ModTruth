"""Machine-readable CLI with explicit scientific and process status separation."""
import argparse
import json
import math
from pathlib import Path
import sys

from . import __version__
from . import b2_core, rmt_core
from .bundle import digest, write_json
from .example import create_example
from .guard import InputProblem, inspect
from .summary import summarize

EXIT = {'PASS': 0, 'REJECT': 1, 'INSUFFICIENT': 2, 'EXECUTION_FAILURE': 3, 'USAGE_ERROR': 64}


class UsageError(Exception):
    pass


class Parser(argparse.ArgumentParser):
    def error(self, message):
        raise UsageError(message)


def audit(bundle, evidence, method='RMT', selection='explicit', threshold=None):
    if method not in ('RMT', 'B2'):
        raise ValueError('UNSUPPORTED_METHOD')
    result = dict(schema_version=1, tool_version=__version__, method=method,
                  operation='SOURCE_BOUND_EXPORT_CHECK', execution_status='COMPLETED',
                  shared_support='Identical schema/hash/composite-key preflight for B2 and RMT',
                  release_state='RESEARCH_SOFTWARE_CANDIDATE',
                  source_authenticity='Caller must independently establish source authenticity; self-hashes are not biological truth')
    try:
        source = inspect(evidence, source=True)
        inspect(bundle, source=False)
    except InputProblem as exc:
        status = 'INSUFFICIENT' if exc.source or exc.code.startswith('UNSUPPORTED_') else 'REJECT'
        result.update(status=status, audit_status=status, raw_status='NOT_RUN_SHARED_PREFLIGHT',
                      raw_result=None, metrics=None,
                      issues=[dict(rule=exc.code, input_role='evidence' if exc.source else 'bundle', detail=exc.detail)])
        return result
    raw = rmt_core.validate(bundle, evidence) if method == 'RMT' else b2_core.check(bundle, evidence)
    normalized = {'ACCEPT': 'PASS', 'REJECT': 'REJECT', 'UNDETERMINED': 'INSUFFICIENT'}[raw['status']]
    metrics = summarize(source, selection, threshold) if normalized == 'PASS' else None
    result.update(status='INSUFFICIENT' if metrics and metrics['metric_status'] == 'INSUFFICIENT' else normalized,
                  audit_status=normalized, raw_status=raw['status'], raw_result=raw, metrics=metrics,
                  metrics_input='Complete verified source ledger after candidate scientific-equivalence PASS' if metrics else None,
                  source_manifest_sha256=digest(Path(evidence) / 'SOURCE_LINKS.json'))
    return result


def main(argv=None):
    parser = Parser(prog='rmt-targeted', description=__doc__)
    parser.add_argument('--version', action='version', version=__version__)
    sub = parser.add_subparsers(dest='command', required=True, parser_class=Parser)
    check = sub.add_parser('check', help='Audit one schema-v1 bundle, then report explicitly conditional metrics')
    check.add_argument('--bundle', required=True)
    check.add_argument('--evidence', required=True)
    check.add_argument('--method', choices=('RMT', 'B2'), default='RMT')
    check.add_argument('--selection', choices=('explicit', 'symmetric', 'one-sided'), default='explicit')
    check.add_argument('--threshold', type=float)
    check.add_argument('--strict-exit', action='store_true')
    check.add_argument('--output', help='Optional new JSON file; existing files are never overwritten')
    example = sub.add_parser('example', help='Create a new wholly fictional bundle/evidence pair')
    example.add_argument('--output', required=True)
    args = None
    try:
        args = parser.parse_args(argv)
        if args.command == 'check':
            if args.selection == 'explicit' and args.threshold is not None:
                raise UsageError('--threshold does not apply to explicit selection')
            if args.selection != 'explicit' and args.threshold is None:
                raise UsageError('--threshold is required for a declared filter')
            if args.threshold is not None and (not math.isfinite(args.threshold) or not (0.5 if args.selection == 'symmetric' else 0) <= args.threshold <= 1):
                raise UsageError('Invalid threshold domain')
        if args.command == 'example':
            result = create_example(args.output)
        else:
            result = audit(args.bundle, args.evidence, args.method, args.selection, args.threshold)
            if args.output:
                write_json(args.output, result)
        print(json.dumps(result, sort_keys=True, allow_nan=False))
        return EXIT[result['status']] if getattr(args, 'strict_exit', False) else 0
    except UsageError as exc:
        result = dict(status='USAGE_ERROR', execution_status='FAILED', error_type=type(exc).__name__, message=str(exc))
    except Exception as exc:
        result = dict(status='EXECUTION_FAILURE', execution_status='FAILED', error_type=type(exc).__name__, message=str(exc))
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return EXIT[result['status']]


if __name__ == '__main__':
    raise SystemExit(main())
