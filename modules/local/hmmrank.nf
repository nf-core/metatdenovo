process HMMRANK {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' } "

    input:

    path hmmrout

    output:

    path "hmmrank"    , emit: x
    path "version.yml", emit: versions

    script:

    """
    #!/usr/bin/env Rscript
    library(optparse)
    library(readr)
    library(dplyr)
    library(dtplyr)
    library(tidyr)
    library(data.table)
    library(stringr)

    # Read all the tblout files
    tlist <- list()
    i <- 1
    for ( tbloutfile in grep('\\.tblout', opt$args, value=TRUE) ) {
      p <- basename(tbloutfile) %>% str_remove('\\..*')
      t <- read_fwf(
        tbloutfile, fwf_cols(content = c(1, NA)), 
        col_types = cols(content = col_character()), 
        comment='#'
      ) %>% 
        filter(! grepl('^ *#', content)) %>%
        separate(
          content, 
          c('accno', 't0', 'profile', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'), 
          '\\s+',  extra='merge', convert = FALSE
        ) %>%
        transmute(accno, profile, evalue = as.double(evalue), score = as.double(score)) %>%
        data.table()
      if ( nrow(t) > 0 ) {
        tlist[[i]] <- t
        i <- i + 1
      }
    }
    if ( length(tlist) > 0  ) {
      tblout <- rbindlist(tlist)
    } else {
      write("No results found, exiting", stderr())
      quit(save = 'no', status = 0)
    }
    setkey(tblout, profile, accno)

    # Filter away entries with lower score than the minimum
    tblout <- lazy_dt(tblout) %>%
        group_by(accno, profile) %>%
        filter(score >= min(min_score, max_score)) %>%
        ungroup() %>%
        select(-max_score) %>%
        as.data.table()

    # write output
    write_tsv(
        tblout %>% arrange(accno, rank),
        'hmmrank.out'
        )
    """
}
