"""Minimal PEP 517 pure-Python wheel builder; no downloaded build dependency.

Only this fixed package layout is supported. This is not a general build system.
"""
import base64
import csv
import hashlib
import io
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parent
DIST = 'rmt_targeted-0.1.0.dist-info'


def get_requires_for_build_wheel(config_settings=None):
    return []


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    metadata = ('Metadata-Version: 2.1\nName: rmt-targeted\nVersion: 0.1.0\n'
                'Summary: Local source-bound export validation and conditional summaries\n'
                'Requires-Python: >=3.11\nDescription-Content-Type: text/markdown\n\n'
                + (ROOT / 'README.md').read_text(encoding='utf-8'))
    payloads = {str(p.relative_to(ROOT / 'src')): p.read_bytes()
                for p in sorted((ROOT / 'src' / 'rmt_targeted').rglob('*.py'))}
    payloads.update({
        DIST + '/METADATA': metadata.encode('utf-8'),
        DIST + '/WHEEL': b'Wheel-Version: 1.0\nGenerator: rmt-targeted-local\nRoot-Is-Purelib: true\nTag: py3-none-any\n',
        DIST + '/entry_points.txt': b'[console_scripts]\nrmt-targeted = rmt_targeted.cli:main\n',
        DIST + '/top_level.txt': b'rmt_targeted\n',
    })
    record = io.StringIO(newline='')
    writer = csv.writer(record, lineterminator='\n')
    for name, data in sorted(payloads.items()):
        value = base64.urlsafe_b64encode(hashlib.sha256(data).digest()).rstrip(b'=').decode('ascii')
        writer.writerow([name, 'sha256=' + value, len(data)])
    writer.writerow([DIST + '/RECORD', '', ''])
    payloads[DIST + '/RECORD'] = record.getvalue().encode('utf-8')
    filename = 'rmt_targeted-0.1.0-py3-none-any.whl'
    target = Path(wheel_directory) / filename
    with target.open('xb') as stream, zipfile.ZipFile(stream, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for name, data in sorted(payloads.items()):
            info = zipfile.ZipInfo(name, date_time=(2020, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o644 << 16
            archive.writestr(info, data)
    return filename
