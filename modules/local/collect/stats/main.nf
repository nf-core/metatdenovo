process COLLECT_STATS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), val(samples), path(trimlogs), path(bblogs), path(idxstats), path(fcs), path(mergetab)

    output:
    path "${meta.id}.overall_stats.tsv.gz", emit: overall_stats
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def read_trimlogs = ""
    if ( trimlogs ) {
        read_trimlogs = """%>%
        mutate(
            d = map(
                sample,
                function(s) {
                    read_tsv(
                        pipe(sprintf("grep 'Reads written (passing filters)' %s_*_trimming_report.txt | sed 's/.*: *//' | sed 's/ .*//' | sed 's/,//g'", s)),
                        col_names = c('n_trimmed'),
                        col_types = 'i'
                    ) %>%
                        mutate(n_trimmed = n_trimmed * 2)
                }
            )
        ) %>%
        unnest(d)
        """
    }

    if ( mergetab ) {
        read_mergetab = """
        mergetab <- read_tsv("${mergetab}", show_col_types = FALSE)
        """
    } else {
        read_mergetab = """
        mergetab <- tibble(sample = character())
        """
    }

    """
    #!/usr/bin/env Rscript

    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(stringr)

    start    <- tibble(sample = c("${samples.join('", "')}"))

    trimming <- tibble(sample = c("${samples.join('", "')}")) ${read_trimlogs}

    idxs <- tibble(fname = Sys.glob('*.idxstats')) %>%
        mutate(
            sample = str_remove(basename(fname), '.idxstats'),
            d = map(
                fname,
                \\(f) {
                    read_tsv(pipe(sprintf("grep -v '^\\\\*' %s", f)), col_names = c('chr', 'length', 'idxs_n_mapped', 'idxs_n_unmapped'), col_types = 'ciii') %>%
                        select(-chr, -length) %>%
                        summarise(idxs_n_mapped = sum(idxs_n_mapped), idxs_n_unmapped = sum(idxs_n_unmapped))
                }
            )
        ) %>%
        unnest(d) %>%
        select(-fname)

    counts <- read_tsv("${fcs}", col_types = 'cciicicid') %>%
        group_by(sample) %>% summarise(n_feature_count = sum(count))


    bbduk <- tibble(sample = character(), n_non_contaminated = integer())
    for ( f in Sys.glob('*.bbduk.log') ) {
        s = str_remove(f, '.bbduk.log')
        bbduk <- bbduk %>%
            union(
                read_tsv(
                    pipe(sprintf("grep 'Result:' %s | sed 's/Result:[ \t]*//; s/ reads.*//' | sed 's/:/\t/'", f)),
                    col_names = c('n_non_contaminated'),
                    col_types = 'i'
                ) %>%
                    mutate(sample = s)
            )
    }
    if ( nrow(bbduk) == 0 ) bbduk <- bbduk %>% select(sample)

    # Add in stats from taxonomy and function
    ${read_mergetab}

    # Write output
    start %>%
        left_join(trimming, by = join_by(sample)) %>%
        left_join(bbduk, by = join_by(sample)) %>%
        left_join(idxs, by = join_by(sample)) %>%
        left_join(counts, by = join_by(sample)) %>%
        left_join(mergetab, by = join_by(sample)) %>%
        arrange(sample) %>%
        write_tsv("${meta.id}.overall_stats.tsv.gz")

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
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.overall_stats.tsv
    gzip ${prefix}.overall_stats.tsv

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
