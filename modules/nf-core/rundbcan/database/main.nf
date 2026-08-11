process RUNDBCAN_DATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.9--pyhdfd78af_0' :
        'quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0' }"

    output:
    path "dbcan_db", emit: dbcan_db
    // No version-reporting output here: storeDir only supports `val`/`path` outputs, not the
    // `tuple`+`eval` shape the official module uses -- and RUNDBCAN_CAZYMEANNOTATION already
    // reports the same dbCAN version, so nothing is lost.

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    run_dbcan database \\
        ${args} \\
        --db_dir dbcan_db \\
        --aws_s3
    """

    stub:
    """
    mkdir -p dbcan_db
    """
}
