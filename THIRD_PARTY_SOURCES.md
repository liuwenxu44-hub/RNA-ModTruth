# Sources, rights and attribution boundary

## Project code and derived summaries

The retained checker, adapter, analysis, example and test code was produced in
the RNA-ModTruth project. No project software license is selected or granted by
this candidate; in particular, there is no MIT license grant. Code-publication
authorization and software licensing are separate decisions. The fictional
example is entirely self-authored and contains zero original study records.

The numerical outputs are project-derived aggregate descriptions of the named
provider releases, not new biological observations. Scientific protocols and source bindings
retain input identities and analysis definitions; they do not grant rights or
claim that omitted files are distributed here.

## Oxford Nanopore processed validation inputs

Provider: Oxford Nanopore Technologies. Source bucket: `ont-open-data`.

- [2025 RNA modification validation data](https://epi2me.nanoporetech.com/rna-mod-validation-data/)
  and prefix `rna-modbase-validation_2025.03`.
- [2026 all-5-mer validation data](https://epi2me.nanoporetech.com/rna-mod-validation-all5mer-2026.07/)
  and prefix `rna-mod-validation-all5mer-2026.07`.
- [AWS Open Data Registry entry](https://registry.opendata.aws/ont-open-data/)
  lists Oxford Nanopore Technologies and CC BY-NC 4.0 for the bucket resource.
- [CC BY-NC 4.0 terms](https://creativecommons.org/licenses/by-nc/4.0/).

The two prefixes' coverage is **inferred from the bucket-level registry entry**,
not a prefix-specific license statement. The reviewed complete prefix listings
had no local LICENSE; the 2026 README did not separately grant redistribution.
This is a provenance/rights limitation, not a claim that public access itself
licenses unrestricted reuse. Noncommercial, attribution and change-indication
conditions must be retained for any applicable reuse; a project software
license must never replace them.

The project changes were processed-file parsing; preservation of sample/run/read
and reference/position identities; ML interval/midpoint representations;
source-linked condition labels; controlled export perturbations; descriptive
aggregation and explicitly post hoc sensitivity. The resulting real read-level
tables and perturbation/recovery bundles are **not distributed here**. A
controlled perturbation of real provider records is not wholly synthetic data.
Exact URLs, original byte counts, hashes, sample associations and source-map
locators remain in `execution/RMT-EXEC-001/canonical_bundles/INPUTS.json`.

## Other tools and specifications

The retained code refers to pysam 0.23.3, samtools and modkit 0.6.4; no third-party
binary, wheel, model or native-tool source is included. Requirements/commands
are provenance, not redistribution of those tools. Consult their own upstream
terms when installing them.

- [SAM optional-field specification](https://samtools.github.io/hts-specs/SAMtags.pdf).
- [modkit 0.6.4 validation documentation](https://github.com/nanoporetech/modkit/blob/v0.6.4/book/src/intro_validate.md).
- [m6Anet upstream repository](https://github.com/GoekeLab/m6anet), inspected
  fixture version `590ec277cb48d61774f0872395099e466022e810`.

Only m6Anet fixture-format/site-count **aggregate inspection results** and the
project-written inspection code are retained. The upstream individual/site
prediction fixture files are excluded. No read-level truth linkage, caller run
or molecule-level validation is claimed for that inspection.
