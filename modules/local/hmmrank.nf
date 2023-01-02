process HMMRANK {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' } "

    input:

    path hmmtargsums

    output:

    path "hmmrank.tsv.gz", emit: hmmrank
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    # Read all the tblout files
    tibble(fname = Sys.glob('*.tbl.gz')) %>%
        mutate(
            profile = basename(fname) %>% str_remove('.tbl.gz'),
            d = purrr::map(
                fname,
                function(f) {
                    read_fwf(f, fwf_cols(content = c(1, NA)), col_types = cols(content = col_character()), comment='#') %>%
                        filter(! grepl('^ *#', content)) %>%
                        separate(
                            content,
                            c('accno', 't0', 'profile_desc', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'),
                            '\\\\s+',  extra='merge', convert = FALSE
                        ) %>%
                        transmute(accno, profile_desc, evalue = as.double(evalue), score = as.double(score))
                }
            )
        ) %>%
        unnest(d) %>%
        select(-fname) %>%
        # Group and calculate a rank based on score and evalue; let ties be resolved by profile in alphabetical order
        group_by(accno) %>%
        arrange(desc(score), evalue, profile) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
        write_tsv('hmmrank.tsv.gz')

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    tidyverse: ", packageVersion('tidyverse'))
        ),
        "versions.yml"
    )
    """
}
