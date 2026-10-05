"""Re-present frozen construct curves; no score, caller, or new curve calculation."""
import csv
from decimal import Decimal
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_REL = 'run_001/results/FILTER_CURVES.tsv'
SOURCE_SHA = '441bdb5accc32aa4dfb2f2d0f5bbed97a03ea28896000e4ed6dee49a944c4411'
THRESHOLDS = tuple(Decimal(x) for x in ('0.5', '0.6', '0.7', '0.8', '0.9', '0.95', '0.99'))
CONSTRUCTS = tuple('DRACH_oligo_' + str(i) for i in range(5))
STATUS = 'REUSED_PRIOR_RESULT'


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    outputs = [ROOT / 'CONSTRUCT_CONSISTENCY.md', ROOT / 'CONSTRUCT_CONSISTENCY_SOURCE.csv']
    if any(path.exists() for path in outputs):
        raise FileExistsError('Existing presentation outputs must not be overwritten')
    parent = ROOT.parent
    source = parent / SOURCE_REL
    inputs = json.loads((ROOT / 'INPUTS.json').read_text())
    if inputs['source_curves_sha256'] != SOURCE_SHA or sha(source) != SOURCE_SHA:
        raise ValueError('SOURCE_NOT_BOUND_TO_FROZEN_PARENT')
    if sha(ROOT / 'PROTOCOL.json') != inputs['protocol_sha256']:
        raise ValueError('SENSITIVITY_PROTOCOL_DRIFT')
    protocol = json.loads((ROOT / 'PROTOCOL.json').read_text())
    if tuple(Decimal(str(x)) for x in protocol['thresholds']) != THRESHOLDS:
        raise ValueError('THRESHOLD_SCOPE_DRIFT')
    indexed = {}
    with source.open(newline='', encoding='utf-8') as stream:
        for number, row in enumerate(csv.DictReader(stream, delimiter='\t'), 2):
            if row['release'] != '2025_DRACH' or row['policy'] != 'symmetric_confidence':
                continue
            identity = row['construct'], Decimal(row['threshold'])
            if identity in indexed:
                raise ValueError('DUPLICATE_SOURCE_SELECTOR')
            indexed[identity] = (number, row)
    expected = {(construct, t) for construct in CONSTRUCTS + ('ALL',) for t in THRESHOLDS}
    if set(indexed) != expected:
        raise ValueError('SOURCE_CONSTRUCT_THRESHOLD_SUPPORT_MISMATCH')
    for number, row in indexed.values():
        if row['state'] != 'ESTIMABLE' or not row['Brier_equal_class']:
            raise ValueError('MISSING_EXISTING_ESTIMATE')
        if int(row['S']) != int(row['S_class_0']) + int(row['S_class_1']):
            raise ValueError('SCORED_DENOMINATOR_RECONCILIATION_FAILURE')
        if int(row['S_class_0']) <= 0 or int(row['S_class_1']) <= 0:
            raise ValueError('MISSING_REQUIRED_CLASS')
        if Decimal(row['threshold']) == Decimal('.5') and row['S'] != row['O']:
            raise ValueError('BASELINE_NOT_EXPLICIT_O')
    for t in THRESHOLDS:
        for field in ('S_class_0', 'S_class_1', 'S', 'N', 'O'):
            if sum(int(indexed[(construct, t)][1][field]) for construct in CONSTRUCTS) != int(indexed[('ALL', t)][1][field]):
                raise ValueError('POOLED_COUNT_RECONCILIATION_FAILURE')

    claims = []
    def source_record(construct, t, metric, table, printed):
        number, row = indexed[(construct, t)]
        selector = dict(release='2025_DRACH', policy='symmetric_confidence', construct=construct, threshold=row['threshold'])
        claims.append(dict(
            claim_id=f'{table}:{construct}:{t}:{metric}', section=table,
            claim_text=f'Existing {metric} at symmetric midpoint threshold {t}, {construct}',
            evidence_status=STATUS, analysis_set='2025_DRACH; all saved explicit O; symmetric midpoint selection',
            metric_or_outcome=metric, estimate=row[metric], printed_value=printed,
            uncertainty='NOT_APPLICABLE_DESCRIPTIVE; ML interval sensitivity is reported separately',
            test_and_correction='NONE; no significance test or independent-replication claim',
            source_file=SOURCE_REL, source_commit_or_version='RMT-METHODS-TARGETED-001 frozen run_001',
            source_sha256_when_frozen=SOURCE_SHA, source_row_number=str(number),
            field_selector=json.dumps(dict(where=selector, field=metric), sort_keys=True, separators=(',', ':')),
            figure_or_table=table, citation_keys='FROZEN_LOCAL_RESULT_ONLY', verification_status='PASS',
            construct=construct, threshold=row['threshold'], truth_unit='provider construct-condition label',
            metric_unit='read-site', score_definition='ML midpoint p=(m+0.5)/256; explicit O',
            N_class_0=row['N_class_0'], O_class_0=row['O_class_0'], S_class_0=row['S_class_0'],
            N_class_1=row['N_class_1'], O_class_1=row['O_class_1'], S_class_1=row['S_class_1'],
        ))
        return printed

    def format_value(construct, threshold, metric, table, places=None):
        value = indexed[(construct, threshold)][1][metric]
        printed = f'{Decimal(value):.{places}f}' if places is not None else f'{int(value):,}'
        return source_record(construct, threshold, metric, table, printed)

    def grid(metric, table, places=None):
        lines = ['| 构建体 | ' + ' | '.join(f'{t:.2f}' for t in THRESHOLDS) + ' |',
                 '| --- | ' + ' | '.join('---:' for _ in THRESHOLDS) + ' |']
        for construct in CONSTRUCTS:
            label = 'Oligo ' + construct.rsplit('_', 1)[1]
            lines.append('| ' + label + ' | ' + ' | '.join(format_value(construct, t, metric, table, places) for t in THRESHOLDS) + ' |')
        return '\n'.join(lines)

    increasing = []
    for construct in CONSTRUCTS:
        a, b, c = (Decimal(indexed[(construct, Decimal(t))][1]['Brier_equal_class']) for t in ('.9', '.95', '.99'))
        if a < b < c:
            increasing.append(construct)
    if len(increasing) != 5:
        raise ValueError('PRIOR_5_OF_5_DIRECTION_CLAIM_NOT_VERIFIED')
    brier_grid = grid('Brier_equal_class', 'T1_equal_class_Brier', 3)
    control_grid = grid('S_class_0', 'T2_retained_control')
    modified_grid = grid('S_class_1', 'T3_retained_modified')
    pooled = ['| 阈值 | 类别等权 Brier | 观测组成 Brier | 保留对照位点 | 保留修饰条件位点 |',
              '| ---: | ---: | ---: | ---: | ---: |']
    for t in THRESHOLDS:
        values = [format_value('ALL', t, metric, 'T4_pooled_reference', places)
                  for metric, places in [('Brier_equal_class', 3), ('Brier_observed', 3), ('S_class_0', None), ('S_class_1', None)]]
        pooled.append('| ' + f'{t:.2f}' + ' | ' + ' | '.join(values) + ' |')
    # Count/direction is a verification of fifteen existing cells, not a new curve.
    checked = [dict(where=dict(release='2025_DRACH', policy='symmetric_confidence', construct=construct,
                              threshold=indexed[(construct, Decimal(t))][1]['threshold']), field='Brier_equal_class')
               for construct in CONSTRUCTS for t in ('.9', '.95', '.99')]
    claim = {k: '' for k in claims[0]}
    claim.update(claim_id='C1_prior_high_threshold_direction', section='Summary',
        claim_text='All 5 of 5 existing construct equal-class Brier values increase from 0.90 to 0.95 and from 0.95 to 0.99',
        evidence_status=STATUS, analysis_set='2025_DRACH; five existing constructs; existing midpoint curves',
        metric_or_outcome='constructs_with_both_adjacent_high_threshold_increases', estimate='5', printed_value='5/5',
        uncertainty='NONE; descriptive directional consistency only', test_and_correction='NONE',
        source_file=SOURCE_REL, source_commit_or_version='RMT-METHODS-TARGETED-001 frozen run_001',
        source_sha256_when_frozen=SOURCE_SHA,
        source_row_number=';'.join(str(indexed[(construct, Decimal(t))][0]) for construct in CONSTRUCTS for t in ('.9', '.95', '.99')),
        field_selector=json.dumps(dict(selectors=checked, verification='count each construct satisfying Brier(.90)<Brier(.95)<Brier(.99); denominator=5'), sort_keys=True, separators=(',', ':')),
        figure_or_table='T1_equal_class_Brier', citation_keys='FROZEN_LOCAL_RESULT_ONLY', verification_status='PASS',
        truth_unit='provider construct-condition label', metric_unit='read-site', score_definition='ML midpoint p=(m+0.5)/256; explicit O')
    claims.append(claim)
    markdown = f'''# 2025 DRACH 构建体曲线：既有结果的完整呈现

JOURNAL_NEUTRAL_ONLY · {STATUS}

已有中点结果中，5/5 个构建体的类别等权 Brier 都在 0.90→0.95 和 0.95→0.99 两步升高。下表同时保留全部七个阈值，避免把这一局部方向误读为全阈值范围内的单调变化。这里呈现的是同一批既有结果，不是新增独立证据；ML 区间和阈值归属的不确定性由本次补充分析的其他表格单独处理。

## 表 1. 各构建体的类别等权 Brier

{brier_grid}

## 表 2. 各构建体保留的对照条件位点数

{control_grid}

## 表 3. 各构建体保留的修饰条件位点数

{modified_grid}

表 1–3 说明：分析单位为读段–位点，标签为提供者声明的构建体条件；它不是逐分子修饰纯度真值。起始集合 O 是已有保存记录中、比对与碱基条件合格且具有显式 ML 概率的位点；按 p=(ML+0.5)/256 取中点，再保留 max(p, 1−p)≥阈值的记录。类别等权 Brier 是保留对照与保留修饰条件两类平均平方误差的各二分之一加权，而不是按两类保留数量加权。阈值 0.50 保留全部 O；后续各列的两类实际计分分母分别见表 2 和表 3。Brier 显示至小数点后三位，计数为整数；源文件中的全精度值保存在随附 CSV 中。未计算置信区间、P 值或多重检验校正；表中方向一致性不等同于独立生物学重复、外部验证或整个人群表现。Oligo 0–4 分别对应原始标识 DRACH_oligo_0–4，未删减构建体或更换编号。

## 汇总参照：全部五个构建体合并

此处 ALL 是合并参照，**不是第六个构建体**，也不参与 5/5 的计数。合并类别等权值先在每个条件类中按其保留位点求均值，再给予两个类相同权重；它不等于对五个构建体的类别等权值取等权平均。观测组成 Brier 则按实际保留位点数加权。

{chr(10).join(pooled)}

## 来源与核验

全部数值直接选自冻结结果 `{SOURCE_REL}`，选择条件为 release=2025_DRACH、policy=symmetric_confidence、全部五个原始构建体及单列的 ALL 参照、全部七个既定阈值。未读取原始 BAM、重算概率、重新抽样或新增 leave-one-construct-out 分析。

- 源结果 SHA-256：`{SOURCE_SHA}`。
- 科学来源绑定 SHA-256：`{sha(ROOT / "INPUTS.json")}`。
- 每个显示值及 5/5 方向陈述的精确字段、行号、选择条件与全精度值见 `CONSTRUCT_CONSISTENCY_SOURCE.csv`。
- 已核对 42 个唯一来源行、各阈值的两类计数与总计数，以及逐构建体高阈值方向；这里只重新呈现已冻结结果。
'''
    # Explicit full-string generation keeps display rounding separate from exact source values.
    with outputs[1].open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(claims[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(claims)
    with outputs[0].open('x', encoding='utf-8') as stream:
        stream.write(markdown)
    with outputs[1].open(newline='', encoding='utf-8') as stream:
        reloaded = list(csv.DictReader(stream))
    if len(reloaded) != 134 or any(row['verification_status'] != 'PASS' for row in reloaded):
        raise ValueError('PRESENTATION_SOURCE_CROSSWALK_FAILED')
    print(json.dumps(dict(status='PASS', evidence_status=STATUS, source_rows=42, crosswalk_rows=len(claims),
                          constructs_increasing_at_both_high_threshold_steps=len(increasing),
                          output_sha256={p.name: sha(p) for p in outputs}), sort_keys=True))


if __name__ == '__main__':
    main()
