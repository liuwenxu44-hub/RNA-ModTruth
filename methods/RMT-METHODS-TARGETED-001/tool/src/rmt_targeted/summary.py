"""Conditional midpoint summaries over a verified, complete saved source ledger.

Selections change the stated metric subset, never certify a truncated full bundle.
"""
from collections import Counter

from .guard import TRUTH_KEY, key


def summarize(data, selection='explicit', threshold=None):
    c, rows, truth, _, _ = data
    if selection not in ('explicit', 'symmetric', 'one-sided'):
        raise ValueError('UNSUPPORTED_SELECTION')
    if selection == 'symmetric' and (threshold is None or not .5 <= threshold <= 1):
        raise ValueError('SYMMETRIC_THRESHOLD_MUST_BE_IN_0_5_TO_1')
    if selection == 'one-sided' and (threshold is None or not 0 <= threshold <= 1):
        raise ValueError('ONE_SIDED_THRESHOLD_MUST_BE_IN_0_TO_1')
    labels = {key(r, TRUTH_KEY): int(r['truth_label']) for r in truth}
    constructs = sorted({r['contig'] for r in rows})
    counts = {('ALL', y): Counter() for y in (0, 1)}
    counts.update({(construct, y): Counter() for construct in constructs for y in (0, 1)})
    for r in rows:
        y = labels[key(r, TRUTH_KEY)]
        eligible = r['alignment_status'] == 'MATCHED_READ_SITE'
        explicit = eligible and bool(r['p_mid'])
        p = float(r['p_mid']) if explicit else None
        keep = explicit and (selection == 'explicit' or (max(p, 1 - p) >= threshold if selection == 'symmetric' else p >= threshold))
        for group in (('ALL', y), (r['contig'], y)):
            g = counts[group]
            g['opportunities'] += 1
            g['N'] += eligible
            g['O'] += explicit
            g['S'] += keep
            if not eligible:
                g['excluded_alignment_or_base_scope'] += 1
            elif not explicit:
                g['eligible_missing_probability'] += 1
            if explicit and not keep:
                g['excluded_by_declared_selection'] += 1
            if keep:
                g['Brier_sum'] += (p - y) ** 2
                g['errors'] += int((p >= .5) != bool(y))
    class_rows = []
    for (construct, y), g in sorted(counts.items()):
        s = g['S']
        class_rows.append(dict(
            construct=construct, truth_class=y,
            **{k: g[k] for k in ('opportunities', 'N', 'O', 'S', 'excluded_alignment_or_base_scope',
                                'eligible_missing_probability', 'excluded_by_declared_selection')},
            S_over_O=s / g['O'] if g['O'] else None,
            S_over_N=s / g['N'] if g['N'] else None,
            Brier=g['Brier_sum'] / s if s else None,
            risk=g['errors'] / s if s else None,
            metric_status='PASS' if s else 'INSUFFICIENT',
        ))
    g0, g1 = counts[('ALL', 0)], counts[('ALL', 1)]
    total = g0 + g1
    s = total['S']
    equal = lambda field: sum(g[field] / g['S'] for g in (g0, g1)) / 2 if g0['S'] and g1['S'] else None
    equal_brier, equal_risk = equal('Brier_sum'), equal('errors')
    fixed = [counts[(construct, y)] for construct in constructs for y in (0, 1)]
    fixed_supported = bool(fixed) and all(g['S'] for g in fixed)
    return dict(
        metric_status='PASS' if s else 'INSUFFICIENT', selection=selection, threshold=threshold,
        conditional_scope='Endpoint-eligible explicit ML-midpoint records retained by the declared policy only',
        population_scope=c['population_claim'],
        full_population_loss_identified=False,
        truth_unit='provider construct-condition label; not molecular purity or native occupancy',
        metric_unit='read-site', decision_threshold=.5,
        score_semantics='Native binary m6A ML interval midpoint; 1-p means non-m6A, not multiclass canonical confidence',
        selection_meaning='Declared midpoint selection, not an inferred historical/native omission mechanism',
        counts={k: total[k] for k in ('opportunities', 'N', 'O', 'S', 'excluded_alignment_or_base_scope',
                                     'eligible_missing_probability', 'excluded_by_declared_selection')},
        denominator_definitions=dict(opportunities='saved truth-overlap rows', N='alignment/base eligible',
                                     O='eligible native explicit probabilities', S='retained O'),
        S_over_O=s / total['O'] if total['O'] else None,
        S_over_N=s / total['N'] if total['N'] else None,
        Brier_observed=total['Brier_sum'] / s if s else None,
        risk_observed=total['errors'] / s if s else None,
        Brier_equal_class=equal_brier, risk_equal_class=equal_risk,
        balanced_accuracy=1 - equal_risk if equal_risk is not None else None,
        fixed_construct_equal_class_Brier=sum(g['Brier_sum'] / g['S'] for g in fixed) / len(fixed) if fixed_supported else None,
        fixed_construct_equal_class_risk=sum(g['errors'] / g['S'] for g in fixed) / len(fixed) if fixed_supported else None,
        fixed_construct_support=constructs,
        missing_required_construct_classes=[dict(construct=k, truth_class=y) for (k, y), g in counts.items() if k != 'ALL' and not g['S']],
        missing_strata_policy='Fixed original construct/class support; no dropping or renormalization of absent strata',
        strata=class_rows,
    )
