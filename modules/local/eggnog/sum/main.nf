process EGGNOG_SUM {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-dtplyr=1.3.1 conda-forge::r-data.table=1.14.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:

    tuple val(meta), path(eggnog)
    path(fcs)

    output:

    tuple val(meta), path("${meta.id}.eggnog_summary.tsv.gz") , emit: eggnog_summary
    path "versions.yml"                                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(readr)
    library(tidyr)
    library(stringr)
    library(tidyverse)

    # call the tables into variables
    eggnog <- read_tsv("${eggnog}", show_col_types = FALSE )

    counts <- list.files(pattern = "*.counts.tsv.gz") %>%
        map_df(~read_tsv(.,  show_col_types  = FALSE)) %>%
        mutate(sample = as.character(sample))

    counts %>%
        inner_join(eggnog, by = 'orf') %>%
        group_by(sample) %>%
        drop_na() %>%
        summarise( value = sum(count), .groups = 'drop') %>%
        add_column(database = "eggnog", field = "n") %>%
        relocate(value, .after = last_col()) %>%
        write_tsv('${meta.id}.eggnog_summary.tsv.gz')

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    dtplyr: ", packageVersion('dtplyr')),
            paste0("    data.table: ", packageVersion('data.table')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    purrr: ", packageVersion('purrr')),
            paste0("    tidyr: ", packageVersion('tidyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml")
    """
}
