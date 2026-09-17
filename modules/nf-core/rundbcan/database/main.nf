process RUNDBCAN_DATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Bundles awscli alongside dbcan: this process' storeDir can point at an s3:// path (via
    // --dbcan_dbpath), and Nextflow's AWS Batch executor stages storeDir output to S3 by
    // shelling out to `aws` from inside the task's own container -- which the plain dbcan image
    // doesn't have, addresses #471.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/dbcan_awscli:6aa905a698f38e78' :
        'community.wave.seqera.io/library/dbcan_awscli:ef0587e28e88d850' }"

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
