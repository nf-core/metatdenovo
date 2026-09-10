process COLLECT_FEATURECOUNTSSUMMARY {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b997e8d619c30e5ea23a08d9fb7e4b0c9b441f3187b64d65ff1c0df5e12bba0/data' :
        'community.wave.seqera.io/library/r-base_r-r.utils_r-dplyr_r-readr_pruned:b59bb1a4cfb1196e' }"

    input:
    tuple val(meta), path(summaries)

    output:
    tuple val(meta), path("*.featureCounts.tsv"), emit: unassigned, optional: true
    path "versions.yml"                         , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env Rscript

    library(readr)
    library(dplyr)
    library(stringr)

    # featureCounts' *.summary file is a wide table: one "Status" row per outcome
    # category, one data column per input BAM. Reshape every sample's summary into a
    # long (status, sample, count) table, then write one narrow file per "Unassigned_*"
    # category -- CUSTOM_COLLECTSTATS derives a table column's name from the text
    # between the first and second dot of each fcs file it's handed, so
    # "\${meta.caller}.<Status>.featureCounts.tsv" gives that category its own column.
    # "Assigned" is skipped: it's already the caller's per-sample total, reported as
    # the caller-named column CUSTOM_COLLECTFEATURECOUNTS's own output produces.
    #
    # CUSTOM_COLLECTSTATS reads every fcs file it's given in one combined read_tsv()
    # call, which requires them all to share the same columns -- so each file here
    # carries the same 9 columns CUSTOM_COLLECTFEATURECOUNTS's own counts.tsv.gz does
    # (orf/chr/start/end/strand/length/sample/count/tpm), with only sample and count
    # actually populated.
    summaries <- bind_rows(lapply(Sys.glob('*.featureCounts.tsv.summary'), function(f) {
        read_tsv(f, col_types = cols(Status = col_character(), .default = col_integer())) %>%
            rename(count = 2) %>%
            mutate(sample = str_remove(basename(f), '\\\\.${meta.caller}\\\\.featureCounts\\\\.tsv\\\\.summary\$')) %>%
            select(status = Status, sample, count)
    }))

    for ( s in unique(summaries\$status[str_starts(summaries\$status, 'Unassigned_')]) ) {
        summaries %>%
            filter(status == s) %>%
            transmute(
                orf = NA_character_, chr = NA_character_, start = NA_integer_, end = NA_integer_,
                strand = NA_character_, length = NA_integer_, sample, count, tpm = NA_real_
            ) %>%
            write_tsv(paste0('${meta.caller}.', s, '.featureCounts.tsv'))
    }

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    stringr: ", packageVersion('stringr'))
        ),
        "versions.yml"
    )
    """

    stub:
    """
    touch ${meta.caller}.Unassigned_NoFeatures.featureCounts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        readr: 2.0.0
        dplyr: 1.0.7
        stringr: 1.5.0
    END_VERSIONS
    """
}
