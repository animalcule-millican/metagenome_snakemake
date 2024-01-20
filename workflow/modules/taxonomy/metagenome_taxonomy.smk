configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/modules/taxonomy/taxonomy_config.yml"
workdir: config["working_directory"]

def create_tmpdir(config):
    import random
    import os
    import pickle
    with open(config["tmpdir_pickle"], 'rb') as f:
        adj, aml = pickle.load(f)
    temp_dir_base = config["root_tmpdir"]
    tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    while os.path.exists(tmpdir):
        tmpdir = os.path.join(temp_dir_base, f"{random.choice(adj)}-{random.choice(aml)}")
    os.makedirs(tmpdir, exist_ok=True)
    return tmpdir


def get_sample_names(input_directory):
    import glob
    import os
    sample_names = [os.path.basename(x).split(".")[0] for x in glob.glob(f"{input_directory}/*.fastq.gz")]
    sample_files = glob.glob(f"{input_directory}/*.fastq.gz")
    return sample_names, sample_files

sample_names, sample_files = get_sample_names(config['sample_directory'])

wildcard_constraints:
    sample_name = "|".join(sample_names),
    sample_file = "|".join(sample_files),
    taxa = "|".join(config["taxa"]),
    database = "|".join(config["databases"]),
    output_directory = config["output_directory"],
    reference_directory = config["reference_directory"],
    kmer = "|".join(config["kmer"])



rule all:
    input:
        expand("{output_directory}/taxonomy/{sample_name}_{taxa}_k{kmer}.csv", output_directory = config["output_directory"], sample_name = sample_names, taxa = config["taxa"], kmer = config["kmer"])


rule sketch_metagenome:
    input:
        "{output_directory}/reads/{sample_name}.ecc.fq.gz"
    output:
        "{output_directory}/sketch/{sample_name}.sig.gz"
    threads: 1
    resources:
        mem_mb=30000
    conda:
        "branchwater"
    shell:
        """
        sourmash sketch dna -p k=21,k=31,k=51,abund,scaled=1000 -p k=21,k=31,k=51,abund,scaled=100 --output {output} --name {wildcards.sample_name} {input}
        """

rule fastgather:
    input:
        sample = "{output_directory}/sketch/{sample_name}.sig.gz",
        ref = expand("{reference_directory}/sketch/index/{{taxa}}_sketch.k{{kmer}}.idx", reference_directory = config["reference_directory"])
    output:
        "{output_directory}/gather/fast/{sample_name}_{taxa}_k{kmer}.fastgather.csv"
    threads: 24
    resources:
        mem_mb = 30000
    conda:
        "branchwater"
    shell:
        """
        sourmash scripts fastgather -o {output} -t 10000 -k {wildcards.kmer} -c {threads} {input.sample} {input.ref}
        """

rule gather:
    input:
        sample = "{output_directory}/sketch/{sample_name}.sig.gz",
        ref = expand("{reference_directory}/sketch/index/{{taxa}}_sketch.k{{kmer}}.idx", reference_directory = config["reference_directory"]),
        gather = "{output_directory}/gather/fast/{sample_name}_{taxa}_k{kmer}.fastgather.csv"
    output:
        "{output_directory}/gather/{sample_name}_{taxa}_k{kmer}.gather.csv"
    threads: 1
    resources:
        mem_mb = 20000
    conda:
        "branchwater"
    shell:
        """
        sourmash gather {input.sample} {input.ref} --picklist {input.gather}:match_name:ident -o {output} -k {wildcards.kmer} --threshold-bp 10000
        """

rule taxonomy:
    input:
        gather = "{output_directory}/gather/{sample_name}_{taxa}_k{kmer}.gather.csv",
        lin = expand("{reference_directory}/lineage/{{taxa}}_reference_taxonomy.db", reference_directory = config["reference_directory"])
    output:
        "{output_directory}/taxonomy/{sample_name}_{taxa}_k{kmer}.csv"
    threads: 1
    resources:
        mem_mb = 20480
    conda:
        "branchwater"
    shell:
        """
        sourmash tax metagenome --gather {input.gather} --taxonomy {input.lin} --keep-full-identifiers > {output}
        """