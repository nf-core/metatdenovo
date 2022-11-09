process EGGNOG_TABLE {
    tag '$meta.id'
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gzip=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11 ':
        'quay.io/biocontainers/gzip:1.11 ' }"

    input:
    tuple val(meta), path(eggnog)

    output:
    tuple val(meta), path("eggnogs.tsv.gz"), emit: eggtab
    path "versions.yml"                    , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"

    """
    grep -v '^##' $eggnog |sed 's/^#//' | sed 's/query/orf/' | gzip -c > eggnogs.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$( echo \$(gzip --version 2>&1)| sed 's/.* gzip//')
    END_VERSIONS
    """
}
