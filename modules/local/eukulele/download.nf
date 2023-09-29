process EUKULELE_DOWNLOAD {
    tag "$db"
    label 'process_long'

    conda "bioconda::eukulele=2.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:2.0.5--pyh723bec7_0' :
        'biocontainers/eukulele:2.0.5--pyh723bec7_0' }"

    input:
    tuple val(db), path(directory)

    output:
    path "versions.yml"                           , emit: version
    tuple val("${db}"), path("${directory}/${db}"), emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    pushd $directory

    EUKulele \\
        download \\
        $args \\
        --database $db

    popd
    touch $db

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukulele: \$(echo \$(EUKulele --version 2>&1) | sed 's/Running EUKulele with command line arguments, as no valid configuration file was provided.//; s/The current EUKulele version is//g')
    END_VERSIONS

    """
}
