process HMMRANK {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::r-tidyverse=2.0.0 conda-forge::r-dtplyr=1.3.1 conda-forge::r-data.table=1.14.8"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(tblouts)

    output:
    tuple val(meta), path("*.hmmrank.tsv.gz"), emit: hmmrank
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    #!/usr/bin/env Rscript
    library(readr)
    library(dplyr)
    library(tidyr)
    library(stringr)

    # Read all the tblout files

    read_fwf(c('${tblouts.join("','")}'), fwf_cols(content = c(1, NA)), col_types = cols(content = col_character()), comment='#', id = 'fname') %>%
        filter(! grepl('^ *#', content)) %>%
        separate(
            content,
            c('accno', 't0', 'profile_desc', 't1', 'evalue', 'score', 'bias', 'f0', 'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7', 'f8', 'f9', 'f10', 'rest'),
            '\\\\s+',  extra='merge', convert = FALSE
        ) %>%
        transmute(profile = basename(fname) %>% str_remove('${prefix}\\\\.') %>% str_remove('.tbl.gz'), accno, profile_desc, evalue = as.double(evalue), score = as.double(score)) %>%
        # Group and calculate a rank based on score and evalue; let ties be resolved by profile in alphabetical order
        group_by(accno) %>%
        arrange(desc(score), evalue, profile) %>%
        mutate(rank = row_number()) %>%
        ungroup() %>%
        write_tsv('${prefix}.hmmrank.tsv.gz')

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
