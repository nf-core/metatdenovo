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
    tuple val("${task.process}"), val('pigz'), eval('pigz --version 2>&1'), emit: versions_pigz, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    gunzip = file.toString() - '.gz'

    """
    unpigz \\
        -c \\
        -p $task.cpus \\
        ${file} > $gunzip
    """

    stub:
    gunzip = file.toString() - '.gz'

    """
    touch ${gunzip}
    """
}
