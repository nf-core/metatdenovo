process SUM_EGGNOG {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
    'quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' } "

    input:

    tuple val(meta), path(eggnog)
    path(fcs)

    output:

    tuple val(meta), path("${meta.id}.eggnog_summary.tsv") , emit: eggnog_summary
    path "versions.yml"                         , emit: versions

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
    eggnog <- read_tsv("eggnogs.tsv.gz", show_col_types = FALSE )

    counts <- list.files(pattern = "*_counts.tsv.gz") %>%
        map_df(~read_tsv(.,  show_col_types  = FALSE))

    counts %>% select(1, 7) %>%
        right_join(eggnog, by = 'orf') %>%
        group_by(sample) %>%
        drop_na() %>%
        count(orf) %>%
        summarise( value = sum(n), .groups = 'drop') %>%
        add_column(database = "eggnog", field = "n_orfs") %>%
        relocate(value, .after = last_col()) %>%
        write_tsv('eggnog_summary.tsv')

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")), paste0("    dplyr: ", packageVersion('dplyr')),
        paste0("    dtplyr: ", packageVersion('dtplyr')), paste0("    data.table: ", packageVersion('data.table')), paste0("    readr: ", packageVersion('readr')),
        paste0("    purrr: ", packageVersion('purrr')), paste0("    tidyr: ", packageVersion('tidyr')), paste0("    stringr: ", packageVersion('stringr')) ),
        "versions.yml")
    """
}
