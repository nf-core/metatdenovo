process SUMTAXONOMY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), val(db), path(taxonomy)
    path feature_counts
    val  taxname

    output:
    tuple val(meta), path("*_summary.tsv.gz") , emit: taxonomy_summary
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    # Read the taxonomy and counts tables
    taxonomy <- read_tsv("${taxonomy}", show_col_types = FALSE )

    counts <- read_tsv("${feature_counts}", show_col_types = FALSE) %>%
        mutate(sample = as.character(sample))

    # Join the two and count the number of ORFs with assigned taxonomy
    counts %>%
        inner_join(taxonomy, by = 'orf') %>%
        group_by(sample) %>%
        summarise(value = sum(count), .groups = 'drop') %>%
        mutate(database = "${db ?: 'userdb'}", field = "${taxname}") %>%
        relocate(value, .after = last_col()) %>%
        write_tsv('${prefix}_summary.tsv.gz')

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    tidyverse: ", packageVersion('tidyverse'))
        ),
        "versions.yml"
    )
    """
}
