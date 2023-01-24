process FORMAT_PRODIGAL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gzip=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11 ':
        'quay.io/biocontainers/gzip:1.11 ' }"

    input:
    tuple val(meta), path (gff)

    output:
    tuple val(meta), path("${gff}.gz"), emit: format_gff
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    sed -i 's/\\(\\(k[0-9]\\+_[0-9]\\+\\).*\\)ID=[0-9]\\+\\(_[0-9]\\+\\)/\\1ID=\\2\\3/g' $gff
    gzip $gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$( echo \$(gzip --version 2>&1)| sed 's/.* gzip//')
    END_VERSIONS
    """
}
