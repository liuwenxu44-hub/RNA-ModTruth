# RNA-ModTruth: minimal code and frozen-result verification candidate

This candidate contains the source-bound export-checking tool, entirely fictional
examples/tests, retained scientific code, and frozen aggregate result tables.
It does **not** contain the native research files or source-derived read-level
inputs. It is not a publicly complete scientific recalculation package, an
independent replication, a caller benchmark, or a submission-ready release.

RNA-ModTruth project code is licensed under the [MIT License](LICENSE).
Third-party data, tools and source materials retain their own terms; this
software license does not grant rights to omitted study records. See
[THIRD_PARTY_SOURCES.md](THIRD_PARTY_SOURCES.md). The repository remains a
code and frozen-aggregate verification resource, not a complete public
read-level scientific recalculation package.

## One portable, offline check

From this directory, with Python 3.11 or later and no extra runtime packages:

```bash
python -B verify_public_candidate.py
```

The command writes nothing. It checks the closed file inventory and SHA-256
identities, then verifies frozen aggregate scope, the eight displayed Table 2
rows against both methods' retained recovery facts, seven working populations,
all symmetric thresholds, and saved sensitivity classifications. It does not
read a BAM, run a caller, resample reads, or recompute a scientific metric.
`SOURCE_FILES.json` identifies scientific/source files, their current hashes,
and their source-version hashes. Public metadata and path-only adaptations are
explicitly marked; scientific result files are byte-unchanged.
`PUBLIC_CANDIDATE_MANIFEST.json` covers all files distributed here.

## Fictional example

The example has invented sequences, labels and identifiers, and no original
study records. Choose an absent output directory outside this source tree:

```bash
PYTHONPATH=methods/RMT-METHODS-TARGETED-001/tool/src python -B -m rmt_targeted example --output /absolute/new/fictional-example
PYTHONPATH=methods/RMT-METHODS-TARGETED-001/tool/src python -B -m rmt_targeted check --bundle /absolute/new/fictional-example/bundle --evidence /absolute/new/fictional-example/evidence --method RMT
```

Use `--method B2` for the other implementation. Optional local wheel-building
instructions and the narrow input contract are in
[the tool README](methods/RMT-METHODS-TARGETED-001/tool/README.md).

## Synthetic test scope and commands

The four suites below use Python's standard library and the retained local
modules. They contain 54 tests: tool 29, methods 6, targeted analysis 9 and
sensitivity 10. These are fictional/toy engineering checks, not replay of the
omitted research records. From this directory:

```bash
RMT_TEST_WORK_ROOT=/absolute/existing/test-work PYTHONPATH=methods/RMT-METHODS-TARGETED-001/tool/src python -B -m unittest discover -s methods/RMT-METHODS-TARGETED-001/tool/tests -p test_tool.py -v
python -B -m unittest discover -s methods/RMT-METHODS-001/code -p test_methods.py -v
python -B -m unittest discover -s methods/RMT-METHODS-TARGETED-001/code -p test_analysis.py -v
python -B -m unittest discover -s methods/RMT-METHODS-TARGETED-001/sensitivity_001/code -p test_intervals.py -v
```

Replace `/absolute/existing/test-work` with an existing writable directory
outside this source tree. The tool suite requires `RMT_TEST_WORK_ROOT` and
creates temporary fictional files there; omitting this variable is an error.

The separate retained `execution/RMT-EXEC-001/tests/test_core.py` suite has
11 synthetic parser tests but is **not** standard-library-only: it imports
`pysam`. The original lock targets Python 3.12 / Linux x86_64 and pins
`pysam==0.23.3`. With that dependency already available, its command is:

```bash
python -B -m unittest discover -s execution/RMT-EXEC-001/tests -p test_core.py -v
```

The four suites contain 54 fictional/toy checks. The read-only verification
command above is dependency-free. The separate 11-test parser suite needs
`pysam` and is not included in that standard-library test count.

## What is and is not reproducible here

- Runnable here: the fictional tool demonstration, the four explicitly scoped
  standard-library suites above, frozen-file verification and aggregate-level
  claim/value checks.
- Supplied: metric/checker/analysis code, scientific protocol definitions,
  input hashes and existing full-precision aggregate outputs. Public metadata
  and file-location guards are adapted for this repository; statistical
  definitions, saved numbers and B2/RMT decisions are unchanged.
- Not supplied: four real saved read-site connection tables (78,414,919 bytes),
  canonical/controlled/recovered record bundles, source read-ID selections,
  native BAM/FASTA/BED inputs, third-party fixture rows, executable binaries,
  dependency wheels, logs, manuscript files, author records and historical
  output generations.
- Therefore not demonstrated by this candidate: complete recalculation of
  first-100 and five-hash-subset results, record-level recovery, confidence
  filtering or ML sensitivity from openly provided research records.

Full scientific replay still needs omitted record bundles, read-level inputs
and replay metadata; supplying the source scripts does not make those inputs
available. The targeted analysis accepts a workspace argument and checks the
exact input hashes. Sensitivity also checks its scientific parent-table hashes. Obtain suitable inputs and
resolve their licensing before attempting a separate full replay. Do not treat
metadata/hash coherence as biological truth, molecule-level purity, natural
occupancy, independent validation, or evidence of superiority over B2.

Original acquisition provenance is retained in
`execution/RMT-EXEC-001/canonical_bundles/INPUTS.json`. This candidate performs
no acquisition and does not supply an automated download command.
