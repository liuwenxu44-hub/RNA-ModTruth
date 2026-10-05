"""Measure actual subprocess wall time and Linux wait4 peak RSS, not estimates."""
import os
import subprocess
import time
from pathlib import Path


def measured(command,log):
    started=time.perf_counter()
    env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    with Path(log).open('x') as stream:
        proc=subprocess.Popen([str(x) for x in command],stdout=stream,stderr=subprocess.STDOUT,env=env)
        _,status,usage=os.wait4(proc.pid,0)
        proc.returncode=os.waitstatus_to_exitcode(status)
    return dict(command=[str(x) for x in command],exit_code=proc.returncode,
        wall_seconds=time.perf_counter()-started,cpu_user_seconds=usage.ru_utime,
        cpu_system_seconds=usage.ru_stime,peak_rss_kib=usage.ru_maxrss,
        peak_rss_definition='Linux wait4 ru_maxrss for process including waited children; KiB',log=str(log))
