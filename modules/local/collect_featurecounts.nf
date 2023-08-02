process COLLECT_FEATURECOUNTS {
    tag "$meta.id"
    label 'process_high'

    conda "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' }"

    input:
    tuple val(meta), path(inputfiles)

    output:
    tuple val(meta), path("*.tsv.gz")   , emit: counts
    path "versions.yml"                 , emit: versions

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

    tibble(f = Sys.glob('*.featureCounts.txt')) %>%
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
        write_tsv("${prefix}_counts.tsv.gz")

        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")), paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    dtplyr: ", packageVersion('dtplyr')), paste0("    data.table: ", packageVersion('data.table')) ), "versions.yml")

        """
}
