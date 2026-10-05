"""Read-only BAM/reference/truth canary. No caller, realignment or score imputation.

MM deltas and CIGAR are decoded here independently and compared to pysam/htslib.
Every overlapping primary-alignment/label opportunity is retained, including deletion,
unknown/implicit probabilities, non-A calls and strand exclusions. Original files are
never edited. Experimental truth is the provider's construct-condition annotation.
"""
import argparse
import collections
import copy
import csv
import gzip
import hashlib
import io
import json
import re
import subprocess
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[1]
MM_HEADER = re.compile(r'^([ACGTUN])([+-])([A-Za-z]+|[0-9]+)([.?]?)$')
FIELDS = ['sample_id', 'run_id', 'read_id', 'contig', 'ref_pos0', 'strand',
          'query_pos0', 'truth_label', 'modification_code', 'score_kind', 'ml_integer',
          'p_low', 'p_high', 'p_mid', 'score_status', 'read_base', 'alignment_status',
          'bam_record_ordinal', 'truth_line', 'selected_for_variants']

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def reference(path):
    sequences = {}
    name = None
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            if name in sequences:
                raise ValueError('DUPLICATE_FASTA_NAME')
            sequences[name] = ''
        elif line.strip():
            if name is None:
                raise ValueError('FASTA_SEQUENCE_WITHOUT_HEADER')
            sequences[name] += line.strip().upper().replace('U', 'T')
    return sequences

def annotations(path, refs):
    labels = collections.defaultdict(list)
    for n, line in enumerate(Path(path).read_text().splitlines(), 1):
        if not line or line.startswith(('#', 'track', 'browser')):
            continue
        fields = line.split('\t')
        if len(fields) < 6:
            raise ValueError(f'BED6_REQUIRED:{n}')
        chrom, start, end, code, _, strand = fields[:6]
        if code not in ('-', 'a', '28871'):
            raise ValueError(f'ANNOTATION_NOT_M6A_BINARY:{n}:{code}')
        start, end = int(start), int(end)
        if chrom not in refs or not 0 <= start < end <= len(refs[chrom]):
            raise ValueError(f'BED_REFERENCE_BOUNDS:{n}')
        if strand not in ('+', '-'):
            raise ValueError(f'EXPLICIT_TRUTH_STRAND_REQUIRED:{n}')
        for pos in range(start, end):
            expected = 'A' if strand == '+' else 'T'
            if refs[chrom][pos] != expected:
                raise ValueError(f'TRUTH_BASE_CONFLICT:{n}:{chrom}:{pos}:{refs[chrom][pos]}')
            labels[chrom].append((pos, strand, int(code != '-'), n))
    for chrom, rows in labels.items():
        keys = [(p, s) for p, s, _, _ in rows]
        if len(keys) != len(set(keys)):
            raise ValueError(f'DUPLICATE_OR_CONFLICTING_TRUTH:{chrom}')
    return labels

def decode_mm(read):
    if not read.has_tag('MM'):
        return {}, {}, 'NO_MM_TAG'
    sequence = read.get_forward_sequence().upper()
    if read.has_tag('MN') and read.get_tag('MN') != len(sequence):
        raise ValueError('STALE_MN_LENGTH')
    probabilities = list(read.get_tag('ML')) if read.has_tag('ML') else None
    index = 0
    output, modes = {}, {}
    for segment in read.get_tag('MM').split(';'):
        if not segment:
            continue
        items = segment.split(',')
        m = MM_HEADER.fullmatch(items[0])
        if not m:
            raise ValueError('BAD_MM_HEADER')
        base, sign, raw_codes, mode = m.groups()
        codes = [int(raw_codes)] if raw_codes.isdigit() else list(raw_codes)
        bases = [i for i, letter in enumerate(sequence) if base == 'N' or letter == base or (base == 'U' and letter == 'T')]
        strand = int(sign == '-') ^ int(read.is_reverse)
        for code in codes:
            if (base, strand, code) in modes:
                raise ValueError('REPEATED_MM_GROUP')
            modes[(base, strand, code)] = mode or 'omitted'
        cursor = -1
        for delta in items[1:]:
            cursor += int(delta) + 1
            if cursor < 0 or cursor >= len(bases) or int(delta) < 0:
                raise ValueError('MM_OFFSET_OUT_OF_RANGE')
            original_pos = bases[cursor]
            qpos = len(sequence)-1-original_pos if read.is_reverse else original_pos
            for code in codes:
                if probabilities is not None and index >= len(probabilities):
                    raise ValueError('ML_TOO_SHORT')
                value = probabilities[index] if probabilities is not None else None
                key = (qpos, base, strand, code)
                if key in output:
                    raise ValueError('DUPLICATE_MM_SITE')
                output[key] = value
                index += 1
    if probabilities is not None and index != len(probabilities):
        raise ValueError('ML_LENGTH_MISMATCH')
    return output, modes, 'MM_PRESENT'

def cigar_pairs(read):
    q, r = 0, read.reference_start
    pairs = []
    for op, length in read.cigartuples or []:
        if op in (0, 7, 8):
            pairs.extend((q+i, r+i) for i in range(length))
            q += length
            r += length
        elif op in (1, 4):
            pairs.extend((q+i, None) for i in range(length))
            q += length
        elif op in (2, 3):
            pairs.extend((None, r+i) for i in range(length))
            r += length
        elif op in (5, 6):
            if op == 6:
                raise ValueError('PADDED_CIGAR_NOT_SUPPORTED')
        else:
            raise ValueError('UNSUPPORTED_CIGAR')
    if q != read.query_length:
        raise ValueError('CIGAR_QUERY_LENGTH_CONFLICT')
    return pairs

def run(sample, bam, fasta, bed, output, role, variant_reads_per_construct):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=False)
    counts = collections.Counter()
    issues = []
    refs = reference(fasta)
    labels = annotations(bed, refs)
    target_positions = {chrom:{r[0] for r in rows} for chrom,rows in labels.items()}
    quick = subprocess.run(['samtools', 'quickcheck', '-v', str(bam)], capture_output=True, text=True)
    if quick.returncode:
        raise ValueError('SAMTOOLS_QUICKCHECK_FAILED:'+quick.stderr)
    reads, mapped, primary_names = set(), set(), set()
    read_site_keys = set()
    selection_counts = collections.Counter()
    selected_names = set()
    record_path = output / 'read_sites.tsv.gz'
    with pysam.AlignmentFile(str(bam), 'rb', threads=2) as source:
        header = source.header.to_dict()
        observed_refs = dict(zip(source.references, source.lengths))
        if observed_refs != {k:len(v) for k,v in refs.items()}:
            raise ValueError('BAM_FASTA_REFERENCE_DICTIONARY_MISMATCH')
        readgroups = {r['ID']: r for r in header.get('RG', [])}
        with record_path.open('xb') as raw:
            with gzip.GzipFile(fileobj=raw, filename='', mode='wb', mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding='utf-8', newline='') as text:
                    writer = csv.DictWriter(text, fieldnames=FIELDS, delimiter='\t', lineterminator='\n')
                    writer.writeheader()
                    for ordinal, read in enumerate(source, 1):
                        counts['alignment_records'] += 1
                        reads.add(read.query_name)
                        if not read.is_unmapped:
                            mapped.add(read.query_name)
                        if read.is_secondary or read.is_supplementary:
                            counts['excluded_nonprimary_records'] += 1
                            continue
                        if read.query_name in primary_names:
                            counts['duplicate_primary_read_ids'] += 1
                        primary_names.add(read.query_name)
                        if read.is_unmapped:
                            counts['excluded_unmapped_primary_records'] += 1
                            continue
                        counts['mapped_primary_records'] += 1
                        counts['reverse_primary_records'] += int(read.is_reverse)
                        rg = read.get_tag('RG') if read.has_tag('RG') else ''
                        if rg not in readgroups:
                            counts['read_group_identity_conflicts'] += 1
                            issues.append([ordinal, 'RG_NOT_IN_HEADER'])
                        runid = next((s[6:] for s in readgroups.get(rg, {}).get('DS', '').split() if s.startswith('runid=')), rg)
                        try:
                            decoded, modes, tag_status = decode_mm(read)
                            # htslib 1.18/1.21 mishandle some empty numeric-code groups.
                            # Verify explicit positions in a temporary projection only.
                            # Original groups, skip semantics, BAM and ML remain intact.
                            projected = read
                            groups = read.get_tag('MM').split(';') if read.has_tag('MM') else []
                            if any(group and ',' not in group for group in groups):
                                counts['empty_MM_group_compatibility_projection_records'] += 1
                                projected = copy.copy(read)
                                projected.set_tag('MM', ''.join(group+';' for group in groups if ',' in group), 'Z')
                            hts = {(q,b,s,c): (None if v < 0 else v)
                                   for (b,s,c), positions in (projected.modified_bases or {}).items()
                                   for q,v in positions}
                            if decoded != hts:
                                counts['MM_decoder_disagreements'] += 1
                                if len(issues) < 100:
                                    issues.append([ordinal, 'MM_DECODER_DISAGREEMENT'])
                            pairs = cigar_pairs(read)
                            if pairs != read.get_aligned_pairs():
                                counts['CIGAR_decoder_disagreements'] += 1
                                if len(issues) < 100:
                                    issues.append([ordinal, 'CIGAR_DECODER_DISAGREEMENT'])
                        except (ValueError, TypeError, IndexError) as error:
                            counts['record_parse_errors'] += 1
                            if len(issues) < 100:
                                issues.append([ordinal, str(error)])
                            continue
                        counts['records_MM_present'] += int(tag_status == 'MM_PRESENT')
                        counts['decoded_m6A_scores_all_query_positions'] += sum(v is not None for (q,b,s,c),v in decoded.items() if b == 'A' and c in ('a',28871))
                        ref_to_query = {r:q for q,r in pairs if r is not None}
                        selected = False
                        if read.query_name in selected_names:
                            selected = True
                        elif selection_counts[read.reference_name] < variant_reads_per_construct:
                            selection_counts[read.reference_name] += 1
                            selected_names.add(read.query_name)
                            selected = True
                        hit = 0
                        for pos, strand, truth, line in labels.get(read.reference_name, []):
                            if pos not in ref_to_query:
                                continue
                            hit += 1
                            counts['truth_overlapping_alignment_sites'] += 1
                            q = ref_to_query[pos]
                            key = (runid, read.query_name, read.reference_name, pos, strand)
                            if key in read_site_keys:
                                counts['duplicate_read_site_keys'] += 1
                            read_site_keys.add(key)
                            aligned_strand = int(strand == '-')
                            actual = read.query_sequence[q].upper() if q is not None else ''
                            row = dict.fromkeys(FIELDS, '')
                            row.update(sample_id=sample, run_id=runid, read_id=read.query_name,
                                contig=read.reference_name, ref_pos0=pos, strand=strand,
                                query_pos0=q if q is not None else '', truth_label=truth,
                                modification_code='a', score_kind='read_site_ML_interval',
                                bam_record_ordinal=ordinal, truth_line=line,
                                read_base=actual, selected_for_variants=int(selected))
                            if q is None:
                                state, alignment = 'alignment_deletion_or_skip', 'NO_QUERY_BASE'
                            elif int(read.is_reverse) != aligned_strand:
                                state, alignment = 'opposite_read_strand', 'STRAND_EXCLUDED'
                            elif actual != ('T' if read.is_reverse else 'A'):
                                state, alignment = 'canonical_base_mismatch', 'NON_A_QUERY_BASE'
                            else:
                                alignment = 'MATCHED_READ_SITE'
                                target = next(((q,'A',aligned_strand,c) for c in ('a',28871) if (q,'A',aligned_strand,c) in decoded), None)
                                if target is not None and decoded[target] is not None:
                                    value = decoded[target]
                                    state = 'explicit'
                                    row.update(ml_integer=value, p_low=value/256,
                                               p_high=(value+1)/256, p_mid=(value+0.5)/256)
                                elif target is not None:
                                    state = 'explicit_position_ML_missing'
                                else:
                                    mode = next((modes[('A',aligned_strand,c)] for c in ('a',28871) if ('A',aligned_strand,c) in modes), None)
                                    state = ('no_tag' if tag_status == 'NO_MM_TAG' else 'm6A_group_absent') if mode is None else ('question_unknown' if mode == '?' else 'dot_low_probability_assumption' if mode == '.' else 'omitted_low_probability_assumption')
                            row.update(score_status=state, alignment_status=alignment)
                            counts['site_'+state] += 1
                            counts['selected_variant_sites'] += int(selected)
                            writer.writerow(row)
                        if hit == 0:
                            counts['mapped_primary_reads_without_truth_overlap'] += 1
                        counts['aligned_query_positions_without_truth'] += sum(q is not None and r is not None and r not in target_positions.get(read.reference_name, set()) for q,r in pairs)
    counts['input_unique_read_ids'] = len(reads)
    counts['mapped_unique_read_ids'] = len(mapped)
    counts['primary_unique_read_ids'] = len(primary_names)
    counts['unique_truth_read_sites'] = len(read_site_keys)
    counts['evaluable_explicit_read_sites'] = counts['site_explicit']
    counts['unevaluable_truth_overlapping_sites'] = counts['truth_overlapping_alignment_sites']-counts['site_explicit']
    for name in ('MM_decoder_disagreements', 'CIGAR_decoder_disagreements', 'record_parse_errors',
                 'read_group_identity_conflicts', 'duplicate_primary_read_ids', 'duplicate_read_site_keys'):
        counts.setdefault(name, 0)
    summary = {'sample_id':sample, 'role':role, 'truth_scope':'provider-declared synthetic construct condition, not native occupancy or independently measured purity',
        'counts':dict(sorted(counts.items())), 'input_files':[{ 'path':str(Path(p).resolve()), 'bytes':Path(p).stat().st_size,'sha256':sha(p)} for p in (bam,fasta,bed)],
        'reference_lengths':observed_refs, 'truth_target_sites':sum(map(len,labels.values())),
        'read_groups':header.get('RG', []), 'provider_programs':header.get('PG', []),
        'variants_selection':'first N primary mapped unique read IDs per construct in original BAM order, before score/label eligibility filtering',
        'variant_reads_per_construct':variant_reads_per_construct, 'selected_read_ids_per_construct':dict(selection_counts),
        'record_output_sha256':sha(record_path), 'record_output':'read_sites.tsv.gz',
        'samtools_quickcheck_exit':quick.returncode, 'pysam_version':pysam.__version__,
        'htslib_version':pysam.__samtools_version__,
        'htslib_verification_scope':'explicit MM/ML positions; empty groups excluded in temporary projection only, original skip-state evidence retained',
        'example_issues':issues, 'status':'PASS' if not issues and not any(counts[k] for k in
            ('duplicate_read_site_keys','duplicate_primary_read_ids','MM_decoder_disagreements',
             'CIGAR_decoder_disagreements','record_parse_errors','read_group_identity_conflicts')) else 'REQUIRES_REVIEW'}
    (output/'summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'sample':sample,'status':summary['status'],'counts':summary['counts']},ensure_ascii=False),flush=True)
    return summary

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    for arg in ('sample','bam','fasta','bed','output','role'):
        p.add_argument('--'+arg, required=True)
    p.add_argument('--variant-reads-per-construct', type=int, default=100)
    a = p.parse_args()
    run(a.sample, a.bam, a.fasta, a.bed, a.output, a.role, a.variant_reads_per_construct)
