configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/modules/read_processing/config.yml"
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
    return sample_names

sample_names = get_sample_names(config['sample_directory'])

wildcard_constraints:
    sample_name = "|".join(sample_names)

rule all:
    input:
        expand("{output_directory}/reads/{sample_name}.ecc.fq.gz", sample_name=config["sample_names"], output_directory=config["output_directory"]),
        expand("{output_directory}/reads/{sample_name}.norm.fq.gz", sample_name=config["sample_names"], output_directory=config["output_directory"])

rule filter_reads:
    input:
        config["sample_directory"] + "/{sample_name}.fastq.gz"
    output:
        "{output_directory}/reads/{sample_name}.filt.fq.gz"
    params:
        rqc = config["rqcfilterdata"],
        tmpdir = create_tmpdir()
    threads: 12
    resources:
        mem_mb = 90000
    conda:
        "bbtools"
    shell:
        """
        rqcfilter2.sh -Xmx90g in={input} out={output} \
        rqcfilterdata={params.rqc} \
        barcodefilter=f \
        removehuman=t \
        removedog=t \
        removecat=t \
        removemouse=t \
        qtrim=r \
        trimq=5 \
        silvalocal=f \
        scafstats={params.tmpdir}/scaffoldStats.txt \
        kmerstats={params.tmpdir}/kmerStats.txt \
        log={params.tmpdir}/status.log \
        filelist={params.tmpdir}/file-list.txt \
        stats={params.tmpdir}/filterStats.txt \
        stats2={params.tmpdir}/filterStats2.txt \
        ihist=null \
        reproduceName={params.tmpdir}/reproduce.sh \
        usetmpdir=t \
        tmpdir={params.tmpdir} \
        merge=f

        rm -rf {params.tmpdir}
        """

rule ecc_reads:
    input:
        "{output_directory}/reads/{sample_name}.filt.fq.gz"
    output:
        "{output_directory}/reads/{sample_name}.ecc.fq.gz"
    params:
        tmpdir = create_tmpdir()
    threads: 16
    resources:
        mem_mb = 30000
    conda:
        "bbtools"
    shell:
        """
        bbcms.sh -Xmx30g in={input} out={output} bits=4 hashes=3 k=31
        """

rule norm_reads:
    input:
        "{output_directory}/reads/{sample_name}.ecc.fq.gz"
    output:
        "{output_directory}/reads/{sample_name}.norm.fq.gz"
    params:
        tmpdir = create_tmpdir()
    threads: 16
    resources:
        mem_mb = 60000
    conda:
        "bbtools"
    shell:
        """
        bbnorm.sh -Xmx60g in={input} out={output} ecc=f usetmpdir=t tmpdir={params.tmpdir}

        rm -rf {params.tmpdir}
        """
