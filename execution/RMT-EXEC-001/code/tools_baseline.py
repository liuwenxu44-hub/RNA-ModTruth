"""B1: official native-file tools plus an explicit, independently written join wrapper.

Native modkit results are executed once per unchanged parent and hash-bound; case
times measure this export-join wrapper, not repeated native modkit validation.
No RMT or B2 decisions imported. Unsupported probability/population semantics are
not official-tool errors. B2 separately supplies the full simple semantic baseline.
"""
import argparse
import collections
import subprocess
from pathlib import Path
import pysam
from bundle import ROOT, load, read_json, digest, write_json


def check(path,evidence):
    c,rows,truth,refs,mapping=load(path)
    ec,source,et,ef,em=load(evidence)
    for n,h in read_json(Path(evidence)/'SOURCE_LINKS.json')['files'].items():
        if digest(Path(evidence)/n)!=h:
            raise ValueError('SOURCE_INTEGRITY')
    issues=[]
    allids=[r['evidence_id'] for r in source]
    def add(rule,ids,loc,evi):
        if ids:
            issues.append(dict(rule=rule,evidence_ids=sorted(set(ids)),locator=loc,evidence_locator=evi))
    native=[]
    for sample in ec['samples']:
        bam=ROOT/sample['source_files'][0]['path']
        result=subprocess.run(['samtools','quickcheck','-v',str(bam)],capture_output=True,text=True)
        with pysam.AlignmentFile(str(bam),'rb') as handle:
            sizes=dict(zip(handle.references,handle.lengths))
            groups=handle.header.to_dict().get('RG',[])
        native.append(dict(bam_sha256=sample['source_files'][0]['sha256'],quickcheck_exit=result.returncode,header_references=sizes,read_groups=len(groups)))
        if result.returncode or sizes!={k:len(v) for k,v in ef.items()}:
            add('NATIVE_REFERENCE',allids,str(bam),'samtools quickcheck/pysam reference dictionary')
    ren=c['contig_aliases']
    canonical=lambda n:ren.get(n,n)
    if len(set(ren.values()))!=len(ren) or {canonical(k):v for k,v in refs.items()}!=ef:
        add('REFERENCE_CONTENT',allids,'reference.fa','official FASTA, content-preserving alias mapping only')
    if c['producer_origin']!=c['consumer_origin']:
        add('COORDINATE_DECLARATION',allids,'contract.json#consumer_origin','explicit producer-origin convention')
    if set(c['truth_join_keys'])!=set(ec['truth_join_keys']):
        add('JOIN_KEY',allids,'contract.json#truth_join_keys','BAM sample identity + BED sample mapping + strand/coordinate')
    if set(c['index_join_keys'])!=set(ec['index_join_keys']) or mapping!=em:
        add('READ_MAP',allids,'read_map.tsv','sample/run-index to native UUID/record locator')
    def key(r,origin,norm):
        return tuple([r['sample_id'],norm(r['contig']),str(int(r['position'])-origin),r['strand'],r['truth_label'],r['source_bed_sha256'],r['source_line']])
    if collections.Counter(key(r,c['truth_origin'],canonical) for r in truth)!=collections.Counter(key(r,0,lambda n:n) for r in et):
        add('BED_SAMPLE_MAPPING',allids,'truth.tsv','official modkit-compatible BED source SHA256/line mapped to sample')
    expected={r['evidence_id']:r for r in source}
    bad_identity=[]
    bad_alignment=[]
    bad_mod=[]
    unknown=[]
    seen=collections.Counter()
    for r in rows:
        eid=r['evidence_id']
        seen[eid]+=1
        if eid not in expected:
            unknown.append(eid)
            continue
        e=expected[eid]
        if any(r[k]!=e[k] for k in ('sample_id','run_id','read_id','read_index')):
            bad_identity.append(eid)
        if (canonical(r['contig']),str(int(r['position'])-c['producer_origin']),r['strand'],r['query_pos0'])!=(e['contig'],e['position'],e['strand'],e['query_pos0']):
            bad_alignment.append(eid)
        if r['modification_code'] not in ('a','28871'):
            bad_mod.append(eid)
    add('READ_IDENTITY',bad_identity,'predictions','BAM UUID/ordinal and verified per-run index mapping')
    add('ALIGNMENT_JOIN',bad_alignment,'predictions','native alignment CIGAR and reference coordinate')
    add('REQUESTED_MOD_CODE',bad_mod,'predictions#modification_code','SAMtags m6A a/28871 rather than inosine 17596')
    add('DUPLICATE_READ_SITE',[eid for eid,n in seen.items() if n>1],'predictions','native source record/site keys must not repeat')
    add('UNLOCATED',unknown,'predictions','source record missing')
    return dict(method='B1',status='UNDETERMINED' if unknown else 'REJECT' if issues else 'ACCEPT',
        issues=issues,rule_counts=dict(collections.Counter(x['rule'] for x in issues)),input_records=len(rows),
        native_checks=native,unsupported=['site-vs-read derived probability semantics','filtered population scope','missing-to-zero score corruption','recovery'],
        scope='official native checks + explicit reference/sample/read-ID/coordinate/mod-code joins; native modkit attestation reported separately')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ('bundle','evidence','output'):
        p.add_argument('--'+k,required=True)
    a=p.parse_args()
    write_json(a.output,check(a.bundle,a.evidence))
