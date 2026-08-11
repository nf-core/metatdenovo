process DBCAN_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(overview)

    output:
    tuple val(meta), path("*.dbcan.tsv.gz"), emit: dbcantsv
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    sed 's/^Gene ID/orf/' ${overview} | gzip -c > ${prefix}.dbcan.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gzip -c /dev/null > ${prefix}.dbcan.tsv.gz
    """
}
