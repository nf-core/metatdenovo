process FORMAT_TAX {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' }"

    input:
    tuple val(meta), path(taxtable)

    output:
    tuple val(meta), path("*.taxonomy_classification.tsv.gz"), emit: tax
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(stringr)

    # create a table with taxonomy categories in each column
    read_tsv(Sys.glob('*.out.gz')) %>%
        select(-1) %>%
        rename(orf = transcript_name) %>%
        group_by(orf) %>%
        filter(max_pid == max(max_pid)) %>%
        ungroup() %>%
        separate(
            full_classification,
            c("domain","phylum", "class", "order", "family", "genus", "species"), 
            sep = ";"
        ) %>%
        mutate(
            domain  = ifelse(is.na(domain)  | domain  == '', 'Uncl.',                     domain),
            phylum  = ifelse(is.na(phylum)  | phylum  == '', sprintf("%s uncl.", domain), phylum),
            class   = ifelse(is.na(class)   | class   == '', sprintf("%s uncl.", phylum), class),
            order   = ifelse(is.na(order)   | order   == '', sprintf("%s uncl.", class),  order),
            family  = ifelse(is.na(family)  | family  == '', sprintf("%s uncl.", order),  family),
            genus   = ifelse(is.na(genus)   | genus   == '', sprintf("%s uncl.", family), genus),
            species = ifelse(is.na(species) | species == '', sprintf("%s uncl.", genus),  species)
        ) %>%
        write_tsv("${prefix}.taxonomy_classification.tsv.gz")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dplyr: ", packageVersion("dplyr")) ), "versions.yml")
    """
}
