process FORMAT_EUKULELE_TAX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

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
    library(tidyr)

    # Create and write a table with taxonomy categories in each column
    read_tsv("${taxtable}") %>%
        select(-1) %>%
        rename(orf = transcript_name) %>%
        group_by(orf) %>%
        filter(max_pid == max(max_pid)) %>%
        ungroup() %>%
        separate(
            full_classification,
            c("domain","phylum", "class", "order", "family", "genus", "species"),
            sep = "\\\\s*;\\\\s*"
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

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion("dplyr"))
        ),
        "versions.yml"
    )
    """
}
