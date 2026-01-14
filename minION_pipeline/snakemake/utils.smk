import os
import pandas as pd
import csv



def read_mapping_file(filename):
    samples = []
    with open(filename, 'r') as file:
        next(file)  # Skip the first line
        for line in file:
            columns = line.strip().split('\t')  # Assuming tab-separated values
            samples.append(columns[0])
    return samples

filename = config["map_file"]
samples = read_mapping_file(filename)
print(samples)

def get_var(bash_var):
    var = os.environ.get(bash_var)
    return int(var)

threads_per_job = get_var('SLURM_CPUS_PER_TASK')
print("threads available:", threads_per_job)

mem = os.environ.get('SLURM_MEM_PER_NODE')
if isinstance(mem, str) and len(mem) > 0:
    total_mem_gb = int(get_var('SLURM_MEM_PER_NODE') / 1024) 
else:
    total_mem_gb = int(get_var('SLURM_MEM_PER_CPU') * threads_per_job  / 1024 )
print("memory available:", total_mem_gb, "Gb")

tmpdir = str(get_var('SLURM_JOB_ID'))
print("jobID:", tmpdir)

wDIR = os.getcwd()
print("working directory:", wDIR)



############ cluster config ###############
def default_clust(cc):
    """adds __default__ values from cluster config to other
       rules lacking those keys"""
    dkeys = set(cc["__default__"].keys())
    rules = [k for k in cc if k != "__default__"]
    for r in rules:
        toadd = dkeys.difference(cc[r].keys())
        if toadd:
            cc[r].update(dict((k, cc["__default__"][k]) for k in toadd))
    return(cc)

######## load clusterFILE  ##########
with open(config["clusterFILE"]) as f:
    clust_conf = default_clust(yaml.safe_load(f))
#########  end of this mysterious section   ##### 

