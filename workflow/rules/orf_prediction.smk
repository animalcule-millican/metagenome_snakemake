import os


configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/config.yml"
workdir: config["working_directory"]

def calc_chunks(fasta):
    import gzip
    import math
    if ".gz" in fasta:
        with gzip.open(fasta, "rt") as f:
            lines = f.readlines()
        return math.ceil((len(lines)/2)/12)
    else:
        with open(fasta, "r") as f:
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

rule gene_pred:
    input:
        "{data_dir}/reads/{sample_name}.fastq.gz"
    output:
        prot = "{data_dir}/proteins/{sample_name}.faa",
        nucl = "{data_dir}/genes/{sample_name}.fna",
        gff = "{data_dir}/orf/{sample_name}.gff"
    params:
        chunks = lambda wildcards: calc_chunks(wildcards.data_dir + "/reads/" + wildcards.sample_name + ".fastq.gz")
    threads: 12
    resources:
        mem_mb = 30000,
        tmpdir = create_tmpdir()
    conda:
        "genepred"
    shell:
        """
        pprodigal -i {input} -o {output.gff} -a {output.prot} -d {output.nucl} -f gff -T 12 -C {params.chunks} -p meta
        rm -rf {resources.tmpdir}
        """

rule mag_genepred:
    input:
        "{data_dir}/bins/{bin_name}.fna.gz"
    output:
        prot = "{data_dir}/proteins/{bin_name}.faa",
        nucl = "{data_dir}/genes/{bin_name}.fna",
        gff = "{data_dir}/orf/{bin_name}.gff"
    params:
        chunks = lambda wildcards: calc_chunks(wildcards.data_dir + "/bins/" + wildcards.bin_name + ".fna.gz")
    threads: 12
    resources:
        mem_mb = 10000,
        tmpdir = create_tmpdir()
    conda:
        "genepred"
    shell:
        """
        pprodigal -i {input} -o {output.gff} -a {output.prot} -d {output.nucl} -f gff -T 12 -C {params.chunks}
        rm -rf {resources.tmpdir}
        """
