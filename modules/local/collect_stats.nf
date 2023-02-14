process COLLECT_STATS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' } "

    input:
    tuple val(meta), val(samples), path(trimlogs), path(bblogs), path(idxstats), path(fcs), path(mergetab)

    output:
    path "${meta.id}_overall_stats.tsv", emit: overall_stats
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    if ( trimlogs ) {
        read_trimlogs = """%>%
        mutate(
            d = map(
                sample,
                function(s) {
                    fread(cmd = sprintf("grep 'Reads written (passing filters)' %s*trimming_report.txt | sed 's/.*: *//' | sed 's/ .*//' | sed 's/,//g'", s)) %>%
                        as_tibble()
                }
            )
        ) %>%
        unnest(d) %>%
        rename(n_trimmed = V1) %>%
        mutate(n_trimmed = n_trimmed*2) %>%
        """
    } else {
        read_trimlogs = "%>%"
    }

    if (mergetab) {
        read_mergetab = """

        mergetab <- list.files(pattern = "*_merged_table.tsv" ) %>%
            map_df(~read_tsv(.,  show_col_types  = FALSE))
        """
    } else {
        read_mergetab = """
        mergetab <- data.frame(sample = character(), stringsAsFactors = FALSE)
        """
    }

    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(dplyr)
    library(readr)
    library(purrr)
    library(tidyr)
    library(stringr)

    TYPE_ORDER = c('n_trimmed', 'n_non_contaminated', 'idxs_n_mapped', 'idxs_n_unmapped', 'n_feature_count')

    # Collect stats for each sample, create a table in long format that can be appended to
    t <- tibble(sample = c("${samples.join('", "')}")) ${read_trimlogs}
    # add samtools idxstats output
        mutate(
            i = map(
                sample,
                function(s) {
                    fread(cmd = sprintf("grep -v '^*' %s*idxstats", s), sep = '\\t', col.names = c('chr', 'length', 'idxs_n_mapped', 'idxs_n_unmapped')) %>%
                    lazy_dt() %>%
                    summarise(idxs_n_mapped = sum(idxs_n_mapped), idxs_n_unmapped = sum(idxs_n_unmapped)) %>%
                    as_tibble()
                }
            )
        ) %>%
        unnest(i) %>%
            pivot_longer(2:ncol(.), names_to = 'm', values_to = 'v') %>%
            union(
            # Total observation after featureCounts
                tibble(file = Sys.glob('*_counts.tsv.gz')) %>%
                mutate(d = map(file, function(f) fread(cmd = sprintf("gunzip -c %s", f), sep = '\\t'))) %>%
                as_tibble() %>%
                unnest(d) %>%
                group_by(sample) %>% summarise(n_feature_count = sum(count), .groups = 'drop') %>%
                pivot_longer(2:ncol(.), names_to = 'm', values_to = 'v')
            )

    # Add in stats from BBDuk, if present
    for ( f in Sys.glob('*.bbduk.log') ) {
        s = str_remove(f, '.bbduk.log')
        t <- t %>% union(
            fread(cmd = sprintf("grep 'Result:' %s | sed 's/Result:[ \\t]*//; s/ reads.*//'", f), col.names = c('v')) %>%
            as_tibble() %>%
            mutate(sample = s, m = 'n_non_contaminated')
        )
    }

    # Add in stats from taxonomy and function
    ${read_mergetab}

    # Write the table in wide format
    t %>%
        mutate(m = parse_factor(m, levels = TYPE_ORDER, ordered = TRUE)) %>%
        arrange(sample, m) %>%
        pivot_wider(names_from = m, values_from = v) %>%
        left_join(mergetab, by = 'sample') %>%
        write_tsv('${prefix}_overall_stats.tsv')

        writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")), paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    dtplyr: ", packageVersion('dtplyr')), paste0("    data.table: ", packageVersion('data.table')), paste0("    readr: ", packageVersion('readr')),
            paste0("    purrr: ", packageVersion('purrr')), paste0("    tidyr: ", packageVersion('tidyr')), paste0("    stringr: ", packageVersion('stringr')) ),
            "versions.yml")
    """
}
