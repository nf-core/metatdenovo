process FORMAT_TAX {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::r-tidyverse=1.3.1 conda-forge::r-data.table=1.14.0 conda-forge::r-dtplyr=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' :
        'quay.io/biocontainers/mulled-v2-508c9bc5e929a77a9708902b1deca248c0c84689:0bb5bee2557136d28549f41d3faa08485e967aa1-0' }"

    input:
    tuple val(meta), path(taxtable)

    output:
    path "*_taxonomy_classification.tsv", emit: tax
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(readr)
    library(dtplyr)
    library(dplyr)
    library(purrr)
    library(tidyr)
    library(stringr)

    setDTthreads($task.cpus)

    # create a table with taxonomy categories in each column
    tax <- list.files(pattern = "*.out") %>%
                map_df(~read.table(.,  sep = "\t", header = TRUE, fill = TRUE))

    tax <- tax[,-c(1,3,5,6,7,8)]

    colnames(tax)[1] <- "orf"

    tax <- tax %>%
            separate('full_classification',c("Domain","Phylum", "Class", "Order", "Family", "Genus", "Species"), ";") %>%
            mutate(
                Domain = ifelse(is.na(Domain)   | Domain == '',  sprintf("%s uncl.", Domain), Domain),
                Phylum = ifelse(is.na(Phylum)   | Phylum == '',  sprintf("%s uncl.", Phylum), Phylum),
                Class = ifelse(is.na(Class)     | Class == '',   sprintf("%s uncl.", str_remove(Phylum, ' unclassified')), Class),
                Order = ifelse(is.na(Order)     | Order == '',   sprintf("%s uncl.", str_remove(Class,  ' unclassified')),  Order),
                Family = ifelse(is.na(Family)   | Family == '',  sprintf("%s uncl.", str_remove(Order,  ' unclassified')),  Family),
                Genus = ifelse(is.na(Genus)     | Genus == '',   sprintf("%s uncl.", str_remove(Family, ' unclassified')), Genus),
                Species = ifelse(is.na(Species) | Species == '', sprintf("%s uncl.", str_remove(Genus, ' unclassified')),  Species)
            ) %>%
            na.omit() %>%
            write_tsv("${meta}_taxonomy_classification.tsv")

    writeLines(c("\\"${task.process}\\":", paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),paste0("    dplyr: ", packageVersion("dplyr")) ), "versions.yml")
    """
}
