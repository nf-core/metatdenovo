process PROKKAGFF2TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.prokka-annotations.tsv.gz"), emit: tsv
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(stringr)

    fread(
        cmd = "zgrep -P '\\t' $gff",
        col.names = c('contig', 'gene_caller', 'feature', 'start', 'end', 'a', 'strand', 'b', 'c')
    ) %>%
        separate_rows(c, sep = ';') %>%
        separate(c, c('k', 'v'), sep = '=') %>%
        pivot_wider(names_from = k, values_from = v) %>%
        select(-a, -b) %>%
        rename(orf = ID) %>%
        rename_all(str_to_lower) %>%
        relocate(sort(colnames(.)[8:ncol(.)]), .after = 7) %>%
        relocate(orf) %>%
        as.data.table() %>%
        write_tsv("${prefix}.prokka-annotations.tsv.gz")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    data.table: ", packageVersion("data.table")),
            paste0("    dtplyr: "    , packageVersion("dtplyr")),
            paste0("    dplyr: "     , packageVersion("dplyr")),
            paste0("    tidyr: "     , packageVersion("tidyr")),
            paste0("    readr: "     , packageVersion("readr"))
        ),
        "versions.yml"
    )
    """
}
