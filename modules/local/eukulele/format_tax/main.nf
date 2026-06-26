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
    path "versions.yml"                                      , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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
        filter(max_score == max(max_score)) %>%
        ungroup() %>%
        separate(
            full_classification,
            c("domain","phylum", "class", "order", "family", "genus", "species"),
            sep = "\\\\s*;\\\\s*"
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
