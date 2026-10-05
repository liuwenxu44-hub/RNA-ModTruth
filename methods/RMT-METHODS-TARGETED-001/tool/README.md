# rmt-targeted 0.1.0

Minimal installable command-line wrapper for the existing schema-v1 source-bound
export checkers, plus explicitly conditional saved-score summaries. Python 3.11+
is required. Runtime and wheel building use only the standard library. No raw
reads, native callers, models, network connection, repository layout or original
research files are required for the fictional example.

This is a narrowly scoped export-checking tool, not a generally applicable
benchmarking framework. RNA-ModTruth project code is licensed under the
[MIT License](../../../LICENSE); third-party data and tools retain their own terms.
When redistributing a wheel built by the retained minimal builder, include the
repository-root LICENSE; the builder does not package that file automatically.
The fictional example contains no original study records. No recovery command
or source-adapter/caller is exposed.

## Install and run

From a local copy of this `tool` directory, using a Python environment with pip:

```bash
python -m pip wheel --no-index --no-deps --no-build-isolation --wheel-dir wheelhouse .
python -m pip install --no-index --no-deps wheelhouse/rmt_targeted-0.1.0-py3-none-any.whl
rmt-targeted example --output /absolute/new/fictional-example
rmt-targeted check --bundle /absolute/new/fictional-example/bundle --evidence /absolute/new/fictional-example/evidence --method RMT
rmt-targeted check --bundle /absolute/new/fictional-example/bundle --evidence /absolute/new/fictional-example/evidence --method B2 --selection symmetric --threshold 0.9 --strict-exit
```

Create `wheelhouse` first. The builder and `--output` never overwrite an existing
file. `example --output` requires a nonexistent directory. `check` prints one JSON
document to stdout; `--output NEW.json` additionally writes the same document.
Installed commands and `python -m rmt_targeted` work from any working directory.

## Narrow input contract

Supply a candidate bundle and a separately authenticated source ledger. Each
contains `contract.json`, the declared prediction TSV/TSV.GZ parts, `truth.tsv`,
`reference.fa`, and `read_map.tsv`. Evidence additionally requires
`SOURCE_LINKS.json` with exact SHA-256 hashes for all and only those declared
files. Self-consistent hashes do **not** establish truth or external provenance.
Authenticating the ledger is the operator's responsibility. The package does not
reopen native BAM/BED files or infer omitted probability semantics.

Only the registered `schema_version=1`,
`endpoint=READ_SITE_PROBABILITY_EVALUATION` format is supported. Source coordinates
must be zero-based; candidates may use the documented equivalent declared
coordinate/bijective naming conventions accepted by both checker copies.
Required prediction fields, truth fields, map fields and contract fields are
defined explicitly in `rmt_targeted.guard`. The shipped example is an executable
small specification of this schema. Extra scientific fields are not validated by
the legacy decisions; this is not an extensible arbitrary-table validator.

Stable keys are tuples, not row order or concatenated labels:

- Read map: `(sample_id, run_id, read_index)`.
- Truth join: `(sample_id, contig, position, strand)`.
- Trusted source locator: `(sample_id, run_id, source_bam_sha256,
  bam_record_ordinal, query_pos0, contig, position, strand)`.
- `evidence_id` remains an opaque attested record identifier, independently unique
  in the source ledger. It is not parsed into assumed biological identifiers.

Identical and conflicting source/truth/map keys are never silently deduplicated.
Repeated source locator keys with different evidence IDs are also insufficient
source evidence. Candidate repeated records and map defects are handled by the
same legacy decision logic for B2 and RMT, after identical shared schema support.
Missing source files, incomplete hashes, ledger conflicts, source hash drift and
unsupported source contracts produce `INSUFFICIENT`, not invented results.

## Outcomes and process exit codes

JSON `audit_status` is the normalized source-bound decision; `raw_status` and
`raw_result` preserve the historical vocabulary and rule details. Legacy
`ACCEPT`, `REJECT`, and `UNDETERMINED` map to `PASS`, `REJECT`, and `INSUFFICIENT`.
Shared-preflight failures have `raw_status=NOT_RUN_SHARED_PREFLIGHT`.
No retained scores gives overall `status=INSUFFICIENT` even when
`audit_status=PASS`; the empty denominator and null metrics are explicit.

| Outcome | Default exit | `--strict-exit` |
| --- | ---: | ---: |
| PASS | 0 | 0 |
| REJECT | 0 | 1 |
| INSUFFICIENT | 0 | 2 |
| EXECUTION_FAILURE | 3 | 3 |
| USAGE_ERROR | 64 | 64 |

Default exit 0 on a completed rejection is intentional, not a historical bug.
Automation must inspect JSON status or request strict exits. Filesystem write
failure or unexpected execution errors are not scientific rejection. Package
usage errors are not confused with insufficient scientific information.

## Conditional metrics and denominators

Metrics are calculated from the verified complete source ledger **only after**
the candidate export passes. This preserves identical metrics for row ordering,
split/merge and equivalent representation changes. A truncated export claiming
the complete population still fails; it is never accepted merely by requesting a
selection policy. `--selection` changes only the labelled metric subset.

The output reports saved truth-overlap opportunities, eligible `N`, native
explicit-probability `O`, retained `S`, `S/O`, `S/N`, and mutually interpretable
ineligible, eligible-missing, and explicitly-filtered counts. It reports all
construct/class strata; Brier and error risk on retained midpoint scores;
observed-composition and equal-class means; balanced accuracy; and fixed original
construct × class means. Missing required classes/strata produce null fixed-support
metrics, never zero or silently renormalized estimates.

Policies: `explicit` retains all O; `symmetric` retains `max(p, 1-p) >= t` for
`0.5 <= t <= 1`; `one-sided` retains `p >= t`. For this binary m6A score, non-m6A
does not mean multiclass canonical confidence. This CLI reports midpoint policy
metrics only, not latent threshold uncertainty or native omission bounds. The
targeted analysis A–C reports supply those additional descriptive analyses.

These are read-site/provider-condition descriptions, not molecule purity,
occupancy, calibration against molecule-level truth, causal effects, native
omission validation or comparisons between callers/models/releases. Bundles must
not combine incompatible releases. Input `truth_scope` text is not a validation
certificate.

## Implementation reuse and verification

`rmt_core.py` and `b2_core.py` retain only validation-related definitions from
the immutable `methods/RMT-METHODS-001/code/checkers_v2` versions, with relative
serialization imports. Recovery and source-generation paths are deliberately
omitted. Both methods receive the same strict shared preflight and neither
imports the other's decision logic. `PARENT_SOURCE_BINDING.json` records hashes
and exact extraction details; historical files and historical outcomes remain
unchanged. These changes provide engineering usability, not distinct empirical
utility or an RMT advantage.

The standard-library test suite is in `tests/test_tool.py`. Set
`RMT_TEST_WORK_ROOT` to an existing writable directory outside this source tree
and run it with `PYTHONPATH=src python -B -m unittest discover -s tests -v`.
These fictional-fixture checks do not replay the omitted study records.
