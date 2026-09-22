process COLLECT_PROTEINCONSOLIDATE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4b/4b997e8d619c30e5ea23a08d9fb7e4b0c9b441f3187b64d65ff1c0df5e12bba0/data' :
        'community.wave.seqera.io/library/r-base_r-r.utils_r-dplyr_r-readr_pruned:b59bb1a4cfb1196e' }"

    input:
    tuple val(meta), path(inputfiles), path(provenance), path(clusters)

    output:
    tuple val(meta), path("*.counts.tsv.gz"), emit: counts
    path "versions.yml"                     , emit: versions, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript

    library(data.table)
    library(dtplyr)
    library(readr)
    library(dplyr)
    library(stringr)

    setDTthreads($task.cpus)

    # MMseqs2 clusters, no header: representative, member. A representative is its own member.
    clusters <- fread('${clusters}', sep = '\\t', header = FALSE, col.names = c('cluster', 'orf'))

    # Dedupe on ID: a duplicate row would fan the join out and inflate counts.
    provenance <- fread(cmd = "zcat '${provenance}'", sep = '\\t') %>% distinct(ID, .keep_all = TRUE)

    counts <- tibble(f = Sys.glob('*.featureCounts.tsv')) %>%
        mutate(
            d = purrr::map(
                f,
                function(file) {
                    fread(file, sep = '\\t', skip = 1) %>%
                        melt(measure.vars = c(ncol(.)), variable.name = 'sample', value.name = 'count') %>%
                        lazy_dt() %>%
                        mutate(sample = str_remove(sample, '.sorted.bam')) %>%
                        rename(orf = Geneid, chr = Chr, start = Start, end = End, strand = Strand, length = Length) %>%
                        as_tibble()
                }
            )
        ) %>%
        tidyr::unnest(d) %>%
        mutate(orf = str_remove(orf, '^cds\\\\.')) %>%
        select(-f)

    # Coordinates are per locus, so summarised apart from counts; order by locus id for reproducibility.
    cluster_attrs <- clusters %>%
        left_join(counts %>% distinct(orf, chr, start, end, strand, length), by = 'orf') %>%
        left_join(provenance, by = c('orf' = 'ID')) %>%
        arrange(cluster, orf) %>%
        group_by(cluster) %>%
        summarise(
            chr     = paste(chr, collapse = ';'),
            start   = paste(start, collapse = ';'),
            end     = paste(end, collapse = ';'),
            strand  = paste(strand, collapse = ';'),
            # Representative length, not a sum: tpm is normalised by it.
            length  = length[orf == cluster][1],
            callers = paste(sort(unique(unlist(strsplit(paste(callers, collapse = ','), ',')))), collapse = ','),
            n_calls = sum(n_calls),
            n_loci  = n(),
            loci    = paste(orf, collapse = ','),
            .groups = 'drop'
        )

    # Summing member counts is safe only while each read is assigned to one feature
    # (--bbmap_ambiguous/--featurecounts_fraction). tpm is recomputed from the sums.
    # Fail on loci missing from the cluster table: dropping them would silently lose reads.
    orphans <- counts %>% filter(count > 0) %>% distinct(orf) %>% anti_join(clusters, by = 'orf')
    if (nrow(orphans) > 0) {
        stop(sprintf(
            '%d locus/loci with counts are absent from the cluster table, e.g. %s',
            nrow(orphans), paste(head(orphans\$orf, 5), collapse = ', ')
        ))
    }

    counts %>%
        filter(count > 0) %>%
        inner_join(clusters, by = 'orf') %>%
        group_by(cluster, sample) %>%
        summarise(count = sum(count), .groups = 'drop') %>%
        left_join(cluster_attrs, by = 'cluster') %>%
        group_by(sample) %>%
        mutate(tpm = round((count/length)/sum(count/length) * 1e6, 6)) %>%  # round so consistent
        ungroup() %>%
        rename(orf = cluster) %>%
        select(orf, chr, start, end, strand, length, sample, count, tpm, callers, n_calls, n_loci, loci) %>%
        arrange(orf, sample) %>%
        write_tsv("${prefix}.counts.tsv.gz")

    writeLines(
        c(
            "\\"${task.process}\\":",
            paste0("    R: ", paste0(R.Version()[c("major","minor")], collapse = ".")),
            paste0("    dplyr: ", packageVersion('dplyr')),
            paste0("    readr: ", packageVersion('readr')),
            paste0("    stringr: ", packageVersion('stringr')),
            paste0("    dtplyr: ", packageVersion('dtplyr')),
            paste0("    data.table: ", packageVersion('data.table'))
        ),
        "versions.yml"
    )
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.counts.tsv
    gzip ${prefix}.counts.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: 4.1.0
        dplyr: 1.0.7
        dtplyr: 1.1.0
        data.table: 1.14.0
    END_VERSIONS
    """
}
