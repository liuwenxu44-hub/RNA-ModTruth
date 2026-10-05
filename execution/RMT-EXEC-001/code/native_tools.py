"""Actually execute pinned official checks, retaining counts/filter policies/logs."""
import argparse
import json
import re
import subprocess
from pathlib import Path
from bundle import ROOT,read_json,write_json,digest
from measure import measured


def run(evidence,output):
    evidence,output=Path(evidence).resolve(),Path(output).resolve()
    output.mkdir(parents=True,exist_ok=False)
    contract=read_json(evidence/'contract.json')
    tool=ROOT/'.cache/tools/dist_modkit_v0.6.4_cd85862/modkit'
    versions={name:subprocess.run(cmd,capture_output=True).stdout.decode(errors='replace').strip()
              for name,cmd in [('modkit',[str(tool),'--version']),('samtools',['samtools','--version'])]}
    runs=[]
    pair=[]
    for sample in contract['samples']:
        bam=ROOT/sample['source_files'][0]['path']
        bed=ROOT/sample['source_files'][2]['path']
        for f in sample['source_files']:
            if digest(ROOT/f['path'])!=f['sha256']:
                raise ValueError('NATIVE_INPUT_HASH_DRIFT')
        name=sample['sample_id']
        checks=output/('check_tags_'+name)
        command=[tool,'modbam','check-tags',bam,'--mapped-only','--threads','2','--out-dir',checks,
                 '--log-filepath',output/(name+'_check_tags_detail.log')]
        measurement=measured(command,output/(name+'_check_tags.log'))
        measurement.update(type='modkit_check_tags',sample_id=name,bam_sha256=sample['source_files'][0]['sha256'])
        runs.append(measurement)
        pair+=['--bam-and-bed',bam,bed]
    for mode in ('default_confidence_filter','no_confidence_filter'):
        command=[tool,'validate']+pair+['--canonical-base','A','--threads','2','--suppress-progress',
            '--out-filepath',output/(mode+'.tsv'),'--log-filepath',output/(mode+'_detail.log')]
        if mode=='no_confidence_filter':
            command+=['--filter-threshold','0']
        measurement=measured(command,output/(mode+'.log'))
        measurement.update(type='modkit_validate',mode=mode)
        runs.append(measurement)
    write_json(output/'execution.json',dict(versions=versions,runs=runs,
        scope='Official multiclass argmax + canonical/mod balancing; default quantile filter and threshold=0 separately. Not fixed binary .5 evaluation; no --ignore/renormalization.',
        input_contract_sha256=digest(evidence/'contract.json')))
    print(json.dumps({'output':str(output),'exit_codes':[r['exit_code'] for r in runs]}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--evidence',required=True)
    p.add_argument('--output',required=True)
    a=p.parse_args()
    run(a.evidence,a.output)
