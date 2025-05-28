process KOFAMSCAN_UNIQUE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(kofamtsv)

    output:
    tuple val(meta), path("*-uniq.tsv.gz"), emit: kofamuniq
    path "versions.yml"                   , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(stringr)

    read_tsv("$kofamtsv", col_types = 'ccdddc') %>%
        # Select first the best scoring hit
        group_by(orf) %>% filter(score == max(score)) %>% ungroup() %>%
        # then the hit with the smallest evalue
        group_by(orf) %>% filter(evalue == min(evalue)) %>% ungroup() %>%
        # and last the first row
        group_by(orf) %>% filter(row_number() == 1) %>% ungroup() %>%
        write_tsv("${prefix}.kofamscan-uniq.tsv.gz")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    purrr: ", packageVersion('purrr')),
            paste0("    tidyr: ", packageVersion('tidyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gzip -c /dev/null > ${prefix}.unique.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        readr: 2.0.0
        purrr: 0.3.4
        tidyr: 1.1.3
        stringr: 1.4.0
    END_VERSIONS
    """
}
