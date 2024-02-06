rule deepbgc_coassembly:
    input: 
        "{output_directory}/contigs/{project_name}.coassembly.fna.gz"
    output:
        "{output_directory}/deepbgc/coassembly/{project_name}.bgc.gbk",
        "{output_directory}/deepbgc/coassembly/{project_name}.bgc.tsv",
        "{output_directory}/deepbgc/coassembly/{project_name}.full.gbk",
        "{output_directory}/deepbgc/coassembly/{project_name}.pfam.tsv",
        "{output_directory}/deepbgc/coassembly/{project_name}.antismash.json"
    threads: 12
    resources:
        mem_mb = 40000
    conda:
        "deepbgc"
    shell:
        """
        deepbgc pipeline -d clusterfinder_retrained -d clusterfinder_original -d clusterfinder_geneborder -d deepbgc --output "{wildcards.output_directory}/deepbgc/coassembly/{wildcards.project_name}" --label clf_ret --label clf_og --label clf_gb --label deep -c product_activity -c product_class {input}
        """

rule antismash_coassembly:
    input: 
        "{output_directory}/contigs/{project_name}.coassembly.fna.gz"
    output:
        "{output_directory}/bgc/antismash/coassembly/{project_name}/{project_name}.gbk",
        "{output_directory}/bgc/antismash/coassembly/{project_name}/{project_name}.json"
    threads: 12
    resources:
        mem_mb = 40000
    conda:
        "antismash"
    shell:
        """
        antismash --skip-zip-file --allow-long-headers --fullhmmer --clusterhmmer --tigrfam --cc-mibig --asf --cb-general --cb-knownclusters --cb-subclusters --pfam2go --rre --output-dir {wildcards.output_directory}/bgc/antismash/coassembly/{wildcards.project_name} --output-basename {wildcards.project_name} --genefinding-tool prodigal -c {threads} {input}
        """


rule deepbgc_assembly:
    input: 
        "{output_directory}/contigs/{sample_name}.fna.gz"
    output:
        "{output_directory}/deepbgc/{sample_name}.bgc.gbk",
        "{output_directory}/deepbgc/{sample_name}.bgc.tsv",
        "{output_directory}/deepbgc/{sample_name}.full.gbk",
        "{output_directory}/deepbgc/{sample_name}.pfam.tsv",
        "{output_directory}/deepbgc/{sample_name}.antismash.json"
    threads: 16
    resources:
        mem_mb = 40000
    conda:
        "deepbgc"
    shell:
        """
        deepbgc pipeline -d clusterfinder_retrained -d clusterfinder_original -d clusterfinder_geneborder -d deepbgc --output "{wildcards.output_directory}/deepbgc/{wildcards.sample_name}" --label clf_ret --label clf_og --label clf_gb --label deep -c product_activity -c product_class {wildcards.output_directory}/genomes/{wildcards.sample_name}.fna.gz
        """

rule antismash_assembly:
    input: 
        "{output_directory}/contigs/{sample_name}.fna.gz"
    output:
        "{output_directory}/bgc/antismash/{sample_name}/{sample_name}.gbk",
        "{output_directory}/bgc/antismash/{sample_name}/{sample_name}.json"
    threads: 16
    resources:
        mem_mb = 20000
    conda:
        "antismash"
    shell:
        """
        antismash --skip-zip-file --allow-long-headers --fullhmmer --clusterhmmer --tigrfam --cc-mibig --asf --cb-general --cb-knownclusters --cb-subclusters --pfam2go --rre --output-dir {wildcards.output_directory}/bgc/antismash/{wildcards.sample_name} --output-basename {wildcards.sample_name} --genefinding-tool prodigal -c {threads} {wildcards.output_directory}/genome/GCA_003725295.1.genomic.fna.gz
        """