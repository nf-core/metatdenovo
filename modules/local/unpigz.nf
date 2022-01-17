process UNPIGZ {
    tag "$file"
    label 'process_low'

    conda (params.enable_conda ? "pigz=2.3.4" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4':
        'quay.io/biocontainers/pigz:2.3.4' }"

    input:
    path file

    output:
    path "$gunzip",      emit: unzipped
    path "versions.yml", emit: versions

    script:

    gunzip = file.toString() - '.gz'

    """
    unpigz \\
        -c \\
        -p $task.cpus \\
        ${file} > $gunzip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$(pigz --version 2>/dev/null | sed 's/pigz //')
    END_VERSIONS
    """
}
