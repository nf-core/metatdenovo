// Safely quote a Groovy value as a single-quoted R string literal -- without this, a value
// containing a quote or backslash could break the generated R syntax.
def rq(v) {
    return "'" + v.toString().replace('\\', '\\\\').replace("'", "\\'") + "'"
}

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

    # CUSTOM_COLLECTSTATS names a table column from the text between a file's first and
    # second dot, so name each output "<caller>.<Status>.featureCounts.tsv". Skip
    # "Assigned": CUSTOM_COLLECTFEATURECOUNTS's own output already reports it, as the
    # caller-named column.
    summaries <- bind_rows(lapply(Sys.glob('*.featureCounts.tsv.summary'), function(f) {
        read_tsv(f, col_types = cols(Status = col_character(), .default = col_integer())) %>%
            rename(count = 2) %>%
            mutate(sample = str_remove(basename(f), paste0('\\\\.', ${rq(meta.caller)}, '\\\\.featureCounts\\\\.tsv\\\\.summary\$'))) %>%
            select(status = Status, sample, count)
    }))
    # str_remove() returns its input unchanged on no match, so a mismatched suffix would
    # otherwise silently leave the full filename as "sample" instead of failing.
    stopifnot(
        "a *.featureCounts.tsv.summary filename didn't match the expected <sample>.<caller> pattern" =
            ! any(str_detect(summaries\$sample, '\\\\.featureCounts\\\\.tsv\\\\.summary\$'))
    )

    # CUSTOM_COLLECTSTATS reads all fcs files together in one call, so every file must share
    # the same columns as CUSTOM_COLLECTFEATURECOUNTS's own counts table, even though only
    # sample/count are populated here.
    for ( s in unique(summaries\$status[str_starts(summaries\$status, 'Unassigned_')]) ) {
        summaries %>%
            filter(status == s) %>%
            transmute(
                orf = NA_character_, chr = NA_character_, start = NA_integer_, end = NA_integer_,
                strand = NA_character_, length = NA_integer_, sample, count, tpm = NA_real_
            ) %>%
            write_tsv(paste0(${rq(meta.caller)}, '.', s, '.featureCounts.tsv'))
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
