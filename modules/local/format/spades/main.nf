process FORMATSPADES {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::gzip=1.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("rnaspades.format_header.transcript.fa.gz") , emit: assembly
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    cp $assembly rnaspades.fa.gz
    gunzip -c rnaspades.fa.gz  | sed 's/>NODE_\\([0-9]*\\).*/>NODE_\\1/g' | gzip -c > rnaspades.format_header.transcript.fa.gz
    """

    stub:

    """
    gzip -c /dev/null > rnaspades.format_header.transcript.fa.gz
    """
}
