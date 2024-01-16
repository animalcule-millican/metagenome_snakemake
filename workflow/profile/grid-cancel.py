#!/usr/bin/env python3

import sys
import htcondor
import subprocess
from os.path import join


def print_and_exit(s):
    print(s)
    exit()


jobID, UUID, clusterID = sys.argv[1].split("_")

#jobDir = "/home/glbrc.org/millican/.condor_jobs/{}_{}".format(jobID, UUID)
#jobLog = join(jobDir, "condor.log")

cmd = ['condor_rm', clusterID]

try:
    subprocess.run(cmd, check=True)
    print_and_exit(f"JobID {jobID}; running on cluster {clusterID}, has been removed.")
except subprocess.CalledProcessError as e:
    print_and_exit("failed: {}".format(e))
    print_and_exit(f"JobID {jobID}; on cluster {clusterID}, has NOT been removed.")
    


#schedd = htcondor.Schedd()

#try:
#    schedd.act(htcondor.JobAction.Remove, f"ClusterId == {clusterID}")
#except OSError as e:
#    print_and_exit("failed: {}".format(e))

#print_and_exit(f"JobID {jobID}; running on cluster {clusterID}, has been removed.")