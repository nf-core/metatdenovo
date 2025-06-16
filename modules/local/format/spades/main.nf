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
    path "versions.yml"                                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    cp $assembly rnaspades.fa.gz
    gunzip -c rnaspades.fa.gz  | sed 's/>NODE_\\([0-9]*\\).*/>NODE_\\1/g' | gzip > rnaspades.format_header.transcript.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: "\$(gzip --version 2>&1 | grep '^gzip' | sed 's/^gzip \\([0-9.]\\+\\).*/\\1/')"
    END_VERSIONS
    """
}
