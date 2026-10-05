#!/usr/bin/env python3
"""Read-only integrity and communication checks; never recalculate science."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys


M = Path("methods/RMT-METHODS-001")
T = Path("methods/RMT-METHODS-TARGETED-001")
S = T / "sensitivity_001"


def sha(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_json(root: Path, relative: Path | str):
    return json.loads((root / relative).read_text(encoding="utf-8"))


def read_tsv(root: Path, relative: Path | str):
    with (root / relative).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def one(rows, **selector):
    selected = [row for row in rows if all(row.get(key) == value for key, value in selector.items())]
    require(len(selected) == 1, f"Selector is not unique: {selector}")
    return selected[0]


def manifests(root: Path) -> dict:
    public = read_json(root, "PUBLIC_CANDIDATE_MANIFEST.json")
    source = read_json(root, "SOURCE_FILES.json")
    require(public["freeze_status"] == "candidate", "Unexpected public freeze state")
    require(bool(source["files"]), "Scientific source inventory is empty")
    expected = {"PUBLIC_CANDIDATE_MANIFEST.json"}
    require(len({row["path"] for row in public["files"]}) == len(public["files"]), "Duplicate inventory path")
    for manifest in (source, public):
        for row in manifest["files"]:
            relative = Path(row["path"])
            require(not relative.is_absolute() and ".." not in relative.parts, "Unsafe manifest path")
            path = root / relative
            require(path.is_file() and not path.is_symlink(), f"Missing/nonregular file: {relative}")
            require(path.stat().st_size == row["bytes"], f"Size mismatch: {relative}")
            require(sha(path) == row["sha256"], f"SHA-256 mismatch: {relative}")
            if manifest is public:
                expected.add(relative.as_posix())
    actual = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and ".git" not in path.relative_to(root).parts
    }
    require(actual == expected, f"Unclosed inventory: extra={sorted(actual-expected)} missing={sorted(expected-actual)}")
    indexed = {row["path"]: row for row in public["files"]}
    require(all(indexed.get(row["path"]) == row for row in source["files"]), "Source-to-public manifest mismatch")
    return dict(source_files=len(source["files"]), public_files=len(actual), bytes=sum((root / name).stat().st_size for name in actual),
                public_manifest_sha256=sha(root / "PUBLIC_CANDIDATE_MANIFEST.json"))


def verify_table2(root: Path):
    facts = read_json(root, M / "results/server_run_001/facts/FACTS.json")
    checks = read_json(root, M / "results/server_run_001/facts/METHOD_CHECKS.json")
    cases = [
        ("01_sample_reference_BED_mismatch", "Sample/reference mismatch"),
        ("02_declared_zero_one_based_coordinate_error", "Coordinate-base error"),
        ("05_read_index_reuse_across_runs", "Read-index reuse"),
        ("06_duplicate_read_site_counting", "Duplicate read-site"),
    ]
    printed = {
        "DEVELOPMENT": [("0.80590366", 2332), ("NO_RESULT", 0), ("0.07380842", 2332), ("0.07357055", 2333)],
        "EVALUATION": [("0.78082267", 30467), ("NO_RESULT", 0), ("0.07407759", 30467), ("0.07407659", 30468)],
    }
    baseline_printed = {"DEVELOPMENT": ("0.07357893", 2332), "EVALUATION": ("0.07407759", 30467)}
    output = []
    for role in ("DEVELOPMENT", "EVALUATION"):
        canonical = one(facts, role=role, case_id="00_canonical")["recalculated_metrics"]
        base = baseline_printed[role]
        require((f'{canonical["Brier"]:.8f}', canonical["counts"]["evaluable_points"]) == base, "Canonical Table 2 display drift")
        for (case_id, family), expected_before in zip(cases, printed[role]):
            before = one(facts, role=role, case_id=case_id)["recalculated_metrics"]
            before_display = "NO_RESULT" if before["Brier"] is None else f'{before["Brier"]:.8f}'
            require((before_display, before["counts"]["evaluable_points"]) == expected_before, f"Before display drift: {role}/{case_id}")
            for method in ("B2", "RMT"):
                recovered = one(facts, role=role, case_id=f"{case_id}_{method}_recovered")["recalculated_metrics"]
                decision = one(checks, role=role, case_id=case_id, method=method)
                require(decision.get("recovery", {}).get("written") is True, "Recovery was not recorded as written")
                require(recovered["normalized_source_join_sha256"] == canonical["normalized_source_join_sha256"], "Recovered canonical join mismatch")
                require(recovered["Brier"] == canonical["Brier"] and recovered["counts"] == canonical["counts"], "Recovered canonical endpoint mismatch")
            output.append(dict(release=role, error_family=family, before_Brier=before_display, recovered_Brier=base[0],
                               scored_before=expected_before[1], scored_recovered=base[1], canonical_join="RESTORED",
                               methods_verified=["B2", "RMT"]))
    return output


def verify_scope(root: Path):
    metrics = read_json(root, M / "results/server_run_001/robustness/METRICS.json")
    populations = {"full", "original_first100", "seed_11", "seed_29", "seed_47", "seed_71", "seed_101"}
    conditions = {"canonical", "site_mean_as_read", "missing_zero", "high_p_filter"}
    require(len(metrics) == 1092, "Robustness row-count drift")
    require({row["subset"] for row in metrics} == populations, "Working-population scope drift")
    require({row["condition"] for row in metrics} == conditions, "Transformation scope drift")
    require(all(row["analysis_class"] == "POSTHOC_ROBUSTNESS_ANALYSIS" for row in metrics), "Evidence-class drift")
    for release, expected_constructs in (("2025_DRACH", 5), ("2026_all5mer", 32)):
        rows = [row for row in metrics if row["release"] == release]
        constructs = {row["construct"] for row in rows} - {"ALL"}
        require(len(constructs) == expected_constructs, "Construct universe drift")
        for subset in populations:
            for condition in conditions:
                selected = [row for row in rows if row["subset"] == subset and row["condition"] == condition]
                require({row["construct"] for row in selected} == constructs | {"ALL"}, "Missing population/condition/construct cell")
    curves = read_tsv(root, T / "run_001/results/FILTER_CURVES.tsv")
    thresholds = {"0.5", "0.6", "0.7", "0.8", "0.9", "0.95", "0.99"}
    require(len(curves) == 351, "Filter-curve row-count drift")
    for release in ("2025_DRACH", "2026_all5mer"):
        pooled = [row for row in curves if row["release"] == release and row["construct"] == "ALL" and row["policy"] == "symmetric_confidence"]
        require({row["threshold"] for row in pooled} == thresholds and len(pooled) == 7, "Symmetric threshold scope drift")
    for threshold, display in (("0.9", "0.04728"), ("0.95", "0.09329"), ("0.99", "0.16302")):
        row = one(curves, release="2025_DRACH", construct="ALL", policy="symmetric_confidence", threshold=threshold)
        require(f'{float(row["Brier_equal_class"]):.5f}' == display, "Targeted display drift")
    expected_rows = {"CLASS_INTERVALS.tsv": 546, "SUMMARY_INTERVALS.tsv": 273, "FIXED_CONSTRUCT_INTERVALS.tsv": 14,
                     "DIRECTION_CERTIFICATION.tsv": 1760, "DECOMPOSITION_ENCLOSURES.tsv": 28, "PARENT_RECONCILIATION.tsv": 2184}
    for filename, count in expected_rows.items():
        rows = read_tsv(root, S / "run_001/results" / filename)
        require(len(rows) == count, f"Sensitivity row-count drift: {filename}")
    directions = read_tsv(root, S / "run_001/results/DIRECTION_CERTIFICATION.tsv")
    for a, b in (("0.9", "0.95"), ("0.95", "0.99")):
        row = one(directions, release="2025_DRACH", construct="ALL", layer="variable", weighting="equal_class", threshold_from=a, threshold_to=b)
        require(row["status"] == "CERTIFIED_INCREASE", "Saved pooled certificate drift")
    for a, b in zip(("0.5", "0.6", "0.7", "0.8", "0.9", "0.95"), ("0.6", "0.7", "0.8", "0.9", "0.95", "0.99")):
        row = one(directions, release="2026_all5mer", construct="ALL", layer="variable", weighting="equal_class", threshold_from=a, threshold_to=b)
        require(row["status"] == "CERTIFIED_DECREASE", "Saved 2026 certificate drift")
    regression = read_json(root, M / "results/server_run_001/facts/READ_MAP_REGRESSION.json")
    require(len(regression) == 28 and all(row["status"] == "PASS" for row in regression), "Saved regression scope drift")
    return dict(robustness_rows=1092, working_populations=sorted(populations), conditions=sorted(conditions),
                filter_curve_rows=351, symmetric_thresholds=sorted(thresholds, key=float), sensitivity_rows=expected_rows,
                saved_read_map_regression_checks=28)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    try:
        root = args.root.resolve()
        inventory = manifests(root)
        table2 = verify_table2(root)
        scope = verify_scope(root)
        result = dict(status="PASS_FROZEN_AGGREGATE_VERIFICATION_ONLY", inventory=inventory, table2=table2, scope=scope,
                      scope_boundary="Read-only aggregate identity/display verification; not a record-level replay or independent scientific replication",
                      scientific_recalculation=False, research_data_acquisition=False, caller_executed=False,
                      software_license_selected=False, public_complete_recalculation_package=False)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except Exception as error:
        print(json.dumps(dict(status="FAIL", error=str(error), scientific_recalculation=False), indent=2), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
