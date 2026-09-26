process RUNDBCAN_DATABASE {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // storeDir can point at an s3:// path; Nextflow's AWS Batch executor stages that copy by
    // shelling out to `aws` inside the task's own container, which the plain dbcan image lacks.
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
