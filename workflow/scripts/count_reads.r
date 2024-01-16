#!/usr/bin/env Rscript
library("tidyverse")
library("argparse")
library("Rsubread")

# Command line argument parser
Arg_Parser = function(){
    require(argparse)
    # Create a parser
    parser <- ArgumentParser()
    # Add arguments
    parser$add_argument(c('-i',"--input"), nargs="+", required=TRUE, help = "The path to the files to be counted. Multiple files can be provided separated by space. Note: files should be in BAM format.")
    parser$add_argument(c('-a', "--annotation"), required=TRUE, help = "The path to the annotation file. Note: file should be in GTF format.")
    parser$add_argument(c('-t', "--tmpdir"), required=TRUE, help = "The directory to store temporary files")
    parser$add_argument(c('-n', "--Ncpus"), required=FALSE, default=1, help = "The number of threads to use")
    # Parse the arguments
    args <- parser$parse_args()
    return(args)
}

# Parsed arguments
args <- Arg_Parser()
files = args$input
annot = args$annotation
tmpdir = args$tmpdir
threads = args$Ncpus

# Count reads function
count_reads = function(files = files, annot = annot, tmpdir = tmpdir, threads = threads){
    require('Rsubread')
    df.fc = featureCounts(files,
        # annotation
        annot.inbuilt = NULL,
        annot.ext = annot,
        isGTFAnnotationFile = TRUE,
        GTF.featureType = "CDS",
        GTF.attrType = "ID",
        GTF.attrType.extra = c("product", "ko", "cog", "pfam", "tigrfam", "superfamily", "cath_funfam", "ec_number"),
        chrAliases = NULL,
        # level of summarization
        useMetaFeatures = TRUE,
        # overlap between reads and features
        allowMultiOverlap = FALSE,
        minOverlap = 1,
        fracOverlap = 0,
        fracOverlapFeature = 0,
        largestOverlap = FALSE,
        nonOverlap = NULL,
        nonOverlapFeature = NULL,
        # Read shift, extension and reduction
        readShiftType = "upstream",
        readShiftSize = 0,
        readExtension5 = 0,
        readExtension3 = 0,
        read2pos = NULL,
        # multi-mapping reads
        countMultiMappingReads = TRUE,
        # fractional counting
        fraction = FALSE,
        # long reads
        isLongRead = FALSE,
        # read filtering
        minMQS = 0,
        splitOnly = FALSE,
        nonSplitOnly = FALSE,
        primaryOnly = FALSE,
        ignoreDup = FALSE,
        # strandness
        strandSpecific = 0,
        # exon-exon junctions
        juncCounts = FALSE,
        genome = NULL,
        # parameters specific to paired end reads
        isPairedEnd = TRUE,
        countReadPairs = TRUE,
        requireBothEndsMapped = FALSE,
        checkFragLength = FALSE,
        minFragLength = 50,
        maxFragLength = 320,
        countChimericFragments = TRUE,
        autosort = TRUE,
        # number of CPU threads
        nthreads = threads,
        # read group
        byReadGroup = FALSE,
        # report assignment result for each read
        reportReads = NULL,
        reportReadsPath = NULL,
        # miscellaneous
        maxMOp = 10,
        tmpDir = tmpdir,
        verbose = FALSE)
    return(df.fc)
}

# Run Count Reads function
df.count = count_reads(files)

# Write output