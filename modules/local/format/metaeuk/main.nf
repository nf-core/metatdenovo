process FORMAT_METAEUK_GFF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${prefix}_format.gff.gz"), emit: format_gff
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix ?: "${meta.id}"
    cat_input = gff =~ /\.gz$/ ? "gunzip -c ${gff}" : "cat ${gff}"

    // A CDS record with no TCS_ID= attribute must fail rather than silently emit an empty-string
    // id ("ID=;") -- that would still parse as a match downstream (FORMAT_GFF2BED's own ID=
    // extraction), merging unrelated loci into one instead of losing this one record loudly.
    // FORMAT_GFF2BED (modules/local/format/gff2bed/main.nf) guards its own match() the same way --
    // keep both in sync if this guard's contract ever changes.
    """
    $cat_input \\
        | awk 'BEGIN{FS=OFS="\\t"}
            \$3=="CDS" {
                match(\$9, /TCS_ID=[^;]+/)
                if (RSTART == 0) {
                    printf "ERROR: CDS record has no TCS_ID= attribute: %s\\n", \$0 > "/dev/stderr"
                    exit 1
                }
                id = substr(\$9, RSTART+7, RLENGTH-7)
                sub(/_CDS_[0-9]+\$/, "", id)
                print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,"ID="id";"\$9
            }' \\
        | gzip -c > ${prefix}_format.gff.gz
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    gzip -c /dev/null > ${prefix}_format.gff.gz
    """
}
