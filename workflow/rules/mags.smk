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

bin_names = get_mags(config["mag_directory"])

rule sketch_mags:
    input:
        expand("{data_dir}/bins/{bin_name}.fna.gz", bin_name=bin_names, data_dir=config["data_directory"])
    output:
        "{data_dir}/sketch/mags.zip"
    params:
        mag = config["mag_directory"],
        csv = "{data_dir}/mag_list.csv"
    threads: 12
    resources:
        mem_mb = 20000
    conda:
        "branchwater"
    shell:
        """
        scripts/sketch_mags.py {params.mag} {params.csv} {output} {threads}
        """

rule sketch_reads:
    input:
        "{data_dir}/reads/{sample_name}.fastq.gz"
    output:
        "{data_dir}/reads/{sample_name}.sig.gz"
    params:
        param = "k=21,k=31,k=51,scaled=1000,abund",
        sample = "{sample_name}"
    threads: 1
    resources:
        mem_mb = 10000
    conda:
        "branchwater"
    shell:
        """
        sourmash sketch dna {input} -o {output} -p {params.param} --name {params.sample}
        """

rule fastgather_mags:
    input:
        sample = "{data_dir}/reads/{sample_name}.sig.gz",
        mag = "{data_dir}/sketch/mags.zip"
    output:
        "{data_dir}/gather/{sample_name}_k{kmer}.mags.csv"
    params:
        kmer = "{kmer}",
        threads = 12
    resources:
        mem_mb = 20000
    conda:
        "branchwater"
    shell:
        """
        sourmash scripts fastgather -o {output} -t 10000 -k {params.kmer} -c {threads} {input.sample} {input.mag}
        """

rule gather_mags:
    input:
        sample = "{data_dir}/reads/{sample_name}.sig.gz",
        mag = "{data_dir}/sketch/mags.zip",
        gather = "{data_dir}/gather/{sample_name}_k{kmer}.mags.csv"
    output:
        "{data_dir}/gather/{sample_name}_k{kmer}.mag_gather.csv"
    params:
        kmer = "{kmer}",
        threads = 12
    resources:
        mem_mb = 20000
    conda:
        "branchwater"
    shell:
        """
        sourmash gather {input.sample} {input.mag} --picklist {input.gather}:match_name:ident -o {output} -k {params.kmer} --threshold-bp 10000
        """