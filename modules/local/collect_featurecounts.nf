// Import generic module functions
//include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

//params.options = [:]
//options        = initOptions(params.options)


process COLLECT_FEATURECOUNTS {
    //tag "$meta.id"
    label 'process_high'
    
    //publishDir "${params.outdir}",
    //    mode: params.publish_dir_mode,
    //    saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' }"

    //if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0"
    //} else {
    //    container "quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0"
    //}

    input:
        path inputfiles

    output:
        path "*.tsv.gz"    , emit: counts
    //   path "versions.yml", emit: versions

    script:
    //def software = getSoftwareName(task.process)
    def args = task.ext.args ?: ''
    //def prefix = task.ext.prefix ?: "${meta.id}"
    
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
                            group_by(sample) %>% mutate(tpm = r/sum(r) * 1e6) %>% ungroup() %>%
                            select(-r) %>%
                            as_tibble()
                        }
                    )
            ) %>%
            tidyr::unnest(d) %>%
            select(-f) %>%
            write_tsv("summary_table.tsv.gz")

        """
}
