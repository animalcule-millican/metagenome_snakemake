#!/usr/bin/env python3
import sys
import htcondor
from os import makedirs
from os.path import join
from datetime import datetime
from snakemake.utils import read_job_properties
jobscript = sys.argv[1]
job_properties = read_job_properties(jobscript)

UUID = datetime.now().strftime("%d%m%Y%H%M%S")
UUID = str(UUID)
UUID = UUID.replace("_", "-")
jobDir = f"/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/log/{job_properties['rule']}/{job_properties['jobid']}_{UUID}"

makedirs(jobDir, exist_ok=True)

sub = htcondor.Submit(
    {
        "executable": "/bin/bash",
        "arguments": jobscript,
        "max_retries": "1",
        "log": join(jobDir, f"{job_properties['rule']}-{str(datetime.now().strftime('%d%m%Y'))}.$(Cluster).log"),
        "output": join(jobDir, f"{job_properties['rule']}-{str(datetime.now().strftime('%d%m%Y'))}.$(Cluster).out"),
        "error": join(jobDir, f"{job_properties['rule']}-{str(datetime.now().strftime('%d%m%Y'))}.$(Cluster).err"),
        "getenv": "True",
        "request_cpus": job_properties['threads'],
        "batch_name": f"metag-snek_{job_properties['rule']}_{job_properties['jobid']}.$(Cluster)"
    }
)

request_memory = job_properties['resources'].get("mem_mb", None)
if request_memory is not None:
    sub['request_memory'] = str(request_memory)

request_disk = job_properties['resources'].get("disk_mb", None)
if request_disk is not None:
    sub['request_disk'] = str(request_disk)

schedd = htcondor.Schedd()
with schedd.transaction() as txn:
    clusterID = sub.queue(txn)

# print jobid for use in Snakemake
print("{}_{}_{}".format(job_properties['jobid'], UUID, clusterID))