process EGGNOG_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(annotations)

    output:
    tuple val(meta), path("*.emapper.tsv.gz"), emit: emappertsv
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    zgrep -v '^##' ${annotations} | \\
        sed 's/^#// ; s/^query/orf/' | \\
        gzip -c > ${prefix}.emapper.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gzip -c /dev/null > ${prefix}.emapper.tsv.gz
    """
}
