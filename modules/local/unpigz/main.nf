process UNPIGZ {
    tag "$file"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4':
        'biocontainers/pigz:2.3.4' }"

    input:
    tuple val(meta), path(file)

    output:
    tuple val(meta), path("$gunzip") , emit: unzipped
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    gunzip = file.toString() - '.gz'

    """
    unpigz \\
        -c \\
        -p $task.cpus \\
        ${file} > $gunzip
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( echo \$(pigz --version 2>&1) | sed 's/pigz //')
    END_VERSIONS
    """

    stub:
    gunzip = file.toString() - '.gz'
    """
    touch ${gunzip}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: 2.3.4
    END_VERSIONS
    """
}
