process PROKKAGFF2TSV {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3e/3ed83e91da7958208822b05560fd65002a14ae670285e69acd75b252b436efe0/data' :
        'community.wave.seqera.io/library/r-base_r-r.utils_r-dplyr_r-readr_pruned:b59bb1a4cfb1196e' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("*.prokka-annotations.tsv.gz"), emit: tsv
    path "versions.yml"                          , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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
