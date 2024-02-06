

rule maxbin:

rule metabat:

rule vamb:
    input:
        contigs = "{output_directory}/contigs/{project_name}.coassembly.fna.gz",
        bam = expand("{output_directory}/bam/{sample_name}.coassembly.sorted.bam", sample_name = sample_names, output_directory = config["output_directory"])
    output:
        "{output_directory}/mags/vamb/log.txt",
        composition.npz
        abundance.npz
        model.pt
        latent.npz
        clusters.tsv

    shell:
        """
        vamb --outdir path/to/outdir --fasta /path/to/catalogue.fna.gz --bamfiles /path/to/bam/*.bam -o C
        """

rule concoct:
    input:
        contigs = "{contigs}"
    output:
        concoct = "{contigs}.concoct.csv"
    conda:
        "binning_tools"
    shell:
        """
        cut_up_fasta.py {input.contigs} -c 10000 -o 0 --merge_last -b {params.tmpdir}/contigs_10K.bed > {params.tmpdir}/contigs_10K.fa
        concoct_coverage_table.py {params.tmpdir}/contigs_10K.bed {input.bam} > {params.tmpdir}/coverage_table.tsv
        concoct --composition_file {params.tmpdir}/contigs_10K.fa --coverage_file {params.tmpdir}/coverage_table.tsv -b {params.output}/
        merge_cutup_clustering.py {params.output}/clustering_gt1000.csv > {params.output}/clustering_merged.csv
        mkdir {params.bins}
        extract_fasta_bins.py {input.contigs} {params.output}/clustering_merged.csv --output_path {params.bins}/

        perl -pe "s/,/\tconcoct./g;" concoct_clustering_gt1000.csv > concoct.contigs2bin.tsv                
        """

rule contig2bin:
    input:
        mag = "{mag}"
    output:
        contigs = "{mag}.contig2bin.tsv"
    conda:
        "binning_tools"
    shell:
    """
    Fasta_to_Contigs2Bin.sh -e {input.mag} > {output.contigs}
    """