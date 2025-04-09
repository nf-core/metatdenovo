process FORMAT_DIAMOND_TAX_TAXDUMP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(taxfile), path(names), path(nodes), val(ranks)

    output:
    tuple val(meta), path("*.taxonomy-taxdump.tsv.gz"), emit: taxonomy
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    ranks_to_consider = "c('domain', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')"
    if ( ranks ) {
        ranks_to_consider = "c(\"${ranks.tokenize(';').join('\",\"')}\")"
    }

    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    names <- read_tsv(
        pipe("grep 'scientific name' $names | sed 's/\\t|\\t/\\t/g' | sed 's/\\t|\$/\\t/' | cut -f 1,2"),
        col_names = c('taxid', 'name'),
        col_types = 'ic'
    )
    nodes <- read_tsv(
        pipe("cat $nodes | sed 's/\\t|\\t/\\t/g' | sed 's/\\t|\$/\\t/' | cut -f 1,3"),
        col_names = c('taxid', 'rank'),
        col_types = 'ic'
    ) %>%
        filter(rank %in% $ranks_to_consider)

    taxa <- names %>%
        inner_join(nodes, by = join_by(taxid)) %>%
        distinct(name, rank)

    taxassign <- read_tsv("${taxfile}", col_names = c("orf", "taxid", "evalue", "taxonomy"), col_types = 'cidc')

    # Make a translation table from taxonomy string to the well-known ranks above
    tr <- taxassign %>%
        filter(!is.na(taxonomy)) %>%
        distinct(taxonomy) %>%
        transmute(taxonomy, name = taxonomy) %>%
        separate_rows(name, sep = ';') %>%
        distinct(taxonomy, name) %>%
        # Some duplicates might occur for combinations of taxonomy and rank. Allow many-to-many and concatenate the names.
        inner_join(taxa, by = join_by(name), relationship = 'many-to-many') %>%
        group_by(taxonomy, rank) %>%
        arrange(name) %>%
        summarise(name = str_c(name, collapse = '/'), .groups = 'drop') %>%
        pivot_wider(names_from = rank, values_from = name)

    # And left-join with the original table
    taxassign %>%
        left_join(tr, by = join_by(taxonomy)) %>%
        write_tsv("${prefix}.taxonomy-taxdump.tsv.gz")

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
