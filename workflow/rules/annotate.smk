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

rule mags_deepbgc:
    input: 
        "{data_dir}/bins/{bin_name}.fna.gz"
    output:
        "{data_dir}/deepbgc/{bin_name}.bgc.gbk",
        "{data_dir}/deepbgc/{bin_name}.bgc.tsv",
        "{data_dir}/deepbgc/{bin_name}.full.gbk",
        "{data_dir}/deepbgc/{bin_name}.pfam.tsv",
        "{data_dir}/deepbgc/{bin_name}.antismash.json"
    conda:
        "deepbgc"
    shell:
        """
        deepbgc pipeline -d clusterfinder_retrained -d clusterfinder_original -d clusterfinder_geneborder -d deepbgc --output {params} --label clf_ret --label clf_og --label clf_gb --label deep -c product_activity -c product_class {input}
        """

rule sample_deepbgc:
    input:
        "{data_dir}/assembly/{sample_name}.fna.gz"
    output:
        "{data_dir}/deepbgc/{sample_name}.bgc.gbk",
        "{data_dir}/deepbgc/{sample_name}.bgc.tsv",
        "{data_dir}/deepbgc/{sample_name}.full.gbk",
        "{data_dir}/deepbgc/{sample_name}.pfam.tsv",
        "{data_dir}/deepbgc/{sample_name}.antismash.json"
    conda:
        "deepbgc"
    shell:
        """
        deepbgc pipeline --prodigal-meta-mode -d clusterfinder_retrained -d clusterfinder_original -d clusterfinder_geneborder -d deepbgc --output $OUT --label clf_ret --label clf_og --label clf_gb --label deep -c product_activity -c product_class {input}
        """

rule mag_gapseq:
    input:
        "{data_dir}/proteins/{bin_name}.faa"
    output:
        rxn = "{data_dir}/metabolic_model/{bin_name}-all-Reactions.tbl",
        trans = "{data_dir}/metabolic_model/{bin_name}-Transporter.tbl",
        path = "{data_dir}/metabolic_model/{bin_name}-all-Pathways.tbl",
        draft = "{data_dir}/metabolic_model/{bin_name}-draft.RDS",
        wgh = "{data_dir}/metabolic_model/{bin_name}-rxnWeights.RDS",
        xgen = "{data_dir}/metabolic_model/{bin_name}-rxnXgenes.RDS",
    params:
        media = config["media"],
        dirpath = "{data_dir}/metabolic_model"
    conda:
        "gapseq"
    shell:
        """
        gapseq find -p all -M prot -i 30 -c 60 {input} -T {resources.tmpdir} -K {threads} -f {params.dirpath}
        gapseq find-transport -M prot -K {threads} -f {params.dirpath} {input}
        gapseq draft -r {output.rxn} -t {output.rxn} -p {output.rxn} -c {input} -f {params.dirpath}
        gapseq fill -m {output.rxn} -c {output.rxn} -g {output.rxn} -n {params.media} -f {params.dirpath}
        """