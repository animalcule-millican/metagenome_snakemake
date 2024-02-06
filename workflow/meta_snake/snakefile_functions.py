def calc_chunks(fasta):
    import gzip
    import math
    with gzip.open(fasta, "rt") as f:
        lines = f.readlines()
    return math.ceil((len(lines)/2)/12)

def get_mags(dir_path):
    import os
    files = os.listdir(dir_path)
    bin_list = [file.replace(".fna.gz", "") for file in files if file.replace(".fna.gz", "")]
    return bin_list

def create_tmpdir():
    import random
    import os
    import pickle
    with open("/home/glbrc.org/millican/repos/metagenome_snakemake/etc/adj-aml.pkl", 'rb') as f:
        adj, aml = pickle.load(f)
    temp_dir_base = "/home/glbrc.org/millican/TMPDIR"    # Replace with the base path for temporary directories
    # Construct the temporary directory path
    tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    # Check if the directory exists, and find a new combination if it does
    while os.path.exists(tmpdir):
        tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    # Once we find a combination that does not already exist
    # Create the temporary directory
    os.makedirs(tmpdir, exist_ok=True)
    return tmpdir