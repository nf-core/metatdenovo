process FORMAT_LOCUS_CONSOLIDATE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(merged_bed)

    output:
    tuple val(meta), path("${prefix}.gff.gz")        , emit: gff
    tuple val(meta), path("${prefix}.provenance.tsv.gz"), emit: provenance
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    // bedtools merge output (see BEDTOOLS_MERGE's ext.args): chr, start (0-based), end, strand,
    // distinct "<caller>:<original ID>" names, count. ID-inheritance rule: a locus with exactly one
    // distinct contributor inherits that contributor's own original ID verbatim, so a single-caller
    // run is byte-identical to that caller's own per-caller GFF; a locus with several distinct
    // contributors gets a fresh coordinate-derived ID instead.
    """
    awk 'BEGIN {
            FS = OFS = "\\t"
            print "ID", "callers", "n_calls" | "gzip -c > ${prefix}.provenance.tsv.gz"
        }
        {
            start = \$2 + 1
            end   = \$3
            strand = \$4
            n = split(\$5, names, ",")

            if (n == 1) {
                sep = index(names[1], ":")
                id  = substr(names[1], sep + 1)
                callers = substr(names[1], 1, sep - 1)
            } else {
                id = "locus_" \$1 "_" start "_" end "_" strand
                delete seen
                callers = ""
                for (i = 1; i <= n; i++) {
                    sep    = index(names[i], ":")
                    caller = substr(names[i], 1, sep - 1)
                    if (!(caller in seen)) {
                        seen[caller] = 1
                        callers = (callers == "" ? caller : callers "," caller)
                    }
                }
            }

            print \$1, "locus_consolidate", "CDS", start, end, ".", strand, ".", "ID="id | "gzip -c > ${prefix}.gff.gz"
            print id, callers, \$6 | "gzip -c > ${prefix}.provenance.tsv.gz"
        }' ${merged_bed}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -c /dev/null > ${prefix}.gff.gz
    gzip -c /dev/null > ${prefix}.provenance.tsv.gz
    """
}
