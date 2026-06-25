process COLLECT_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b997e8d619c30e5ea23a08d9fb7e4b0c9b441f3187b64d65ff1c0df5e12bba0/data' :
        'community.wave.seqera.io/library/r-base_r-r.utils_r-dplyr_r-readr_pruned:b59bb1a4cfb1196e' }"

    input:
    tuple val(meta), path(inputfiles)

    output:
    tuple val(meta), path("*.counts.tsv.gz"), emit: counts
    path "versions.yml"                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
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
                        rename(orf = Geneid, chr = Chr, start = Start, end = End, strand = Strand, length = Length) %>%
                        group_by(sample) %>%
                        mutate(tpm = r/sum(r) * 1e6) %>% ungroup() %>%
                        select(-r) %>%
                        as_tibble()
                }
            )
        ) %>%
        tidyr::unnest(d) %>%
        # Transdecoder appends "cds." to ORF IDs in the gff file, but does not in the fasta file. Remove to make compatible between tables.
        mutate(orf = str_remove(orf, '^cds\\\\.')) %>%
        select(-f) %>%
        write_tsv("${prefix}.counts.tsv.gz")

        writeLines(
            c(
                "\\"${task.process}\\":",
                paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
                paste0("    dplyr: ", packageVersion('dplyr')),
                paste0("    readr: ", packageVersion('readr')),
                paste0("    stringr: ", packageVersion('stringr')),
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
