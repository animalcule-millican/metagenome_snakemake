configfile: "/home/glbrc.org/millican/repos/metagenome_snakemake/workflow/config/process_reads_config.yaml"
workdir: config["working_directory"]

def create_tmpdir():
    import random
    import os
    import pickle
    with open("/home/glbrc.org/millican/repos/metagenome_snakemake/etc/adj-aml.pkl", 'rb') as f:
        adj, aml = pickle.load(f)
    temp_dir_base = "/home/glbrc.org/millican/TMPDIR"
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
        rqc = config["rqcfilterdata"]
    threads: 12
    resources:
        mem_mb = 90000,
        tmpdir = create_tmpdir()
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
        scafstats={resources.tmpdir}/scaffoldStats.txt \
        kmerstats={resources.tmpdir}/kmerStats.txt \
        log={resources.tmpdir}/status.log \
        filelist={resources.tmpdir}/file-list.txt \
        stats={resources.tmpdir}/filterStats.txt \
        stats2={resources.tmpdir}/filterStats2.txt \
        ihist=null \
        reproduceName={resources.tmpdir}/reproduce.sh \
        usetmpdir=t \
        tmpdir={resources.tmpdir} \
        merge=f

        rm -rf {resources.tmpdir}
        """

rule ecc_reads:
    input:
        "{output_directory}/reads/{sample_name}.filt.fq.gz"
    output:
        "{output_directory}/reads/{sample_name}.ecc.fq.gz"
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
    threads: 16
    resources:
        mem_mb = 60000,
        tmpdir = create_tmpdir()
    conda:
        "bbtools"
    shell:
        """
        bbnorm.sh -Xmx60g in={input} out={output} ecc=f usetmpdir=t tmpdir={resources.tmpdir}

        rm -rf {resources.tmpdir}
        """
