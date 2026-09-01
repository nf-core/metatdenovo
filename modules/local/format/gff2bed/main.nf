process FORMAT_GFF2BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${prefix}.bed"), emit: bed
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix ?: "${meta.id}"
    cat_input = gff =~ /\.gz$/ ? "gunzip -c ${gff}" : "cat ${gff}"

    // BED is 0-based half-open, GFF is 1-based inclusive: start shifts by one, end is unchanged.
    // The "(^|;)ID=" anchor (rather than a bare "ID=") avoids matching inside MetaEuk's
    // Target_ID=/TCS_ID= attributes, which also contain the substring "ID=".
    // The "cds." prefix strip mirrors TIDYVERSE_STRIPCDSPREFIX's existing TransDecoder-ID normalization
    // (modules/local/tidyverse/stripcdsprefix/main.nf), so a locus this caller is the sole contributor to
    // inherits an ID matching that caller's own per-caller counts table exactly.
    """
    $cat_input \\
        | awk -v caller="${meta.caller}" 'BEGIN{FS="\\t"; OFS="\\t"}
            \$3=="CDS" {
                match(\$9, /(^|;)ID=[^;]+/)
                id = substr(\$9, RSTART, RLENGTH)
                sub(/^;?ID=/, "", id)
                sub(/^cds\\./, "", id)
                print \$1, \$4-1, \$5, caller":"id, 0, \$7
            }' \\
        > ${prefix}.bed
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    """
}
