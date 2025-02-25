process FORMAT_DIAMOND_TAX_RANKLIST {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' :
        'biocontainers/mulled-v2-b2ec1fea5791d428eebb8c8ea7409c350d31dada:a447f6b7a6afde38352b24c30ae9cd6e39df95c4-1' }"

    input:
    tuple val(meta), path(taxfile), val(ranks)

    output:
    tuple val(meta), path("*.taxonomy.tsv.gz"), emit: taxonomy
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // Call separate if we have a list of ranks
    def sep = ""
    if ( ranks ) {
        sep = "separate(taxonomy, c(\"${ranks.tokenize(';').join('\",\"')}\"), remove = FALSE, extra = \"merge\", fill = \"right\", sep = \";\") %>%"
    }
    """
    #!/usr/bin/env Rscript

    library(tidyverse)

    read_tsv("${taxfile}", col_names = c("orf", "taxid", "evalue", "taxonomy")) %>%
        ${sep}
        write_tsv("${prefix}.taxonomy.tsv.gz")

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
