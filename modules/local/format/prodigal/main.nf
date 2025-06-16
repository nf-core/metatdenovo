process FORMAT_PRODIGAL_GFF {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path (gff)

    output:
    tuple val(meta), path("${prefix}_format.gff.gz"), emit: format_gff
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    cat_input = gff =~ /\.gz$/ ? "gunzip -c ${gff}" : "cat ${gff}"

    """
    $cat_input | sed 's/^\\([^\\t]\\+\\)\\(.*ID=\\)[0-9]\\+\\(_[0-9]\\+\\)/\\1\\2\\1\\3/' | gzip -c > ${prefix}_format.gff.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: "\$(gzip --version 2>&1 | grep '^gzip' | sed 's/^gzip \\([0-9.]\\+\\).*/\\1/')"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_format.gff
    gzip ${prefix}_format.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: 1.11
    END_VERSIONS
    """
}
