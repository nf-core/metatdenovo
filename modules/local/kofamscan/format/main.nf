process KOFAMSCAN_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(kofamscan_tsv)

    output:
    tuple val(meta), path("*.kofamscan.tsv.gz"), emit: kofamtsv
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "orf	ko	thrshld	score	evalue	ko_definition" | gzip -c > ${prefix}.kofamscan.tsv.gz

    grep -v '#' ${kofamscan_tsv} | cut -f 2-7 | sed 's/\\t"/\\t/' | sed 's/"\$//' | gzip -c >> ${prefix}.kofamscan.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "orf	ko	thrshld	score	evalue	ko_definition" | gzip -c > ${prefix}.kofamscan.tsv.gz
    """
}
