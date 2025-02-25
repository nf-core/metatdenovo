process COLLECT_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(inputfiles)

    output:
    tuple val(meta), path("*.counts.tsv.gz"), emit: counts
    path "versions.yml"                     , emit: versions

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

    tibble(f = Sys.glob('*.featureCounts.tsv')) %>%
        mutate(
            d = purrr::map(
                f,
                function(file) {
                    fread(file, sep = '\\t', skip = 1) %>%
                        melt(measure.vars = c(ncol(.)), variable.name = 'sample', value.name = 'count') %>%
                        lazy_dt() %>%
                        filter(count > 0) %>%
                        mutate(
                            sample = str_remove(sample, '.sorted.bam'),
                            r = count/Length
                        ) %>%
                        rename( orf = Geneid, chr = Chr, start = Start, end = End, strand = Strand, length = Length ) %>%
                        group_by(sample) %>%
                        mutate(tpm = r/sum(r) * 1e6) %>% ungroup() %>%
                        select(-r) %>%
                        as_tibble()
                }
            )
        ) %>%
        tidyr::unnest(d) %>%
        select(-f) %>%
        write_tsv("${prefix}.counts.tsv.gz")

        writeLines(
            c(
                "\\"${task.process}\\":",
                paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                paste0("    dplyr: ", packageVersion('dplyr')),
                paste0("    dtplyr: ", packageVersion('dtplyr')),
                paste0("    data.table: ", packageVersion('data.table'))
            ),
            "versions.yml"
        )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv
    gzip ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        dtplyr: 1.1.0
        data.table: 1.14.0
    END_VERSIONS
    """
}
