process FORMAT_METAEUKFAA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("${prefix}_format.faa.gz"), emit: format_faa
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix    = task.ext.prefix ?: "${meta.id}"
    cat_input = faa =~ /\.gz$/ ? "gunzip -c ${faa}" : "cat ${faa}"

    // MetaEuk fasta headers (target|contig|strand|bitscore|evalue|nExons|lowCoord|...) must match
    // the gff TCS_ID (fields 1,2,3,7) for id joins to work. Fail on headers with too few fields.
    """
    $cat_input \\
        | awk 'BEGIN{FS="|"; OFS="|"}
            /^>/ {
                if (NF < 7) {
                    print "ERROR: MetaEuk fasta header has " NF " fields, expected at least 7: " \$0 > "/dev/stderr"
                    exit 1
                }
                print \$1, \$2, \$3, \$7
                next
            }
            { print }' \\
        | gzip -c > ${prefix}_format.faa.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -c /dev/null > ${prefix}_format.faa.gz
    """
}
