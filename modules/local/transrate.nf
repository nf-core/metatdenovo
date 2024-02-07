process TRANSRATE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::transrate=1.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/transrate:1.0.3--hec16e2b_4':
        'biocontainers/transrate:1.0.3--hec16e2b_4' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("*assemblies_mqc.csv") , emit: assembly_qc
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    transrate \\
        --threads $task.cpus \\
        --assembly $assembly \\
        --output ${prefix}_transrate \\
        $args

    mv ${prefix}_transrate/assemblies.csv ${prefix}_assemblies_mqc.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transrate: \$(transrate --version)
    END_VERSIONS
    """
}
