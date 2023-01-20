process FORMAT_PRODIGAL_OUTPUT {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' }"

    input:
    tuple val(meta), path(inputfile)

    output:
    tuple val(meta), path("*.tsv.gz")   , emit: counts
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(readr)
    library(dplyr)
    library(stringr)

    setDTthreads($task.cpus)
    db_counts <- fread("${prefix}_counts.tsv.gz")

    db_counts$orf <- str_replace_all(db_counts$orf, ".*_", "_")
    db_counts %>%
        unite(orf, chr, orf, sep = '') %>%
        write_tsv("${prefix}_classification.tsv")
