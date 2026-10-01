process METAEUK_DOWNLOAD {
    tag "$db_name"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaeuk:6.a5d39d9--pl5321hf1761c0_2':
        'quay.io/biocontainers/metaeuk:6.a5d39d9--pl5321hf1761c0_2' }"

    input:
    val db_name

    // storeDir allows only val/path outputs, so no version here; METAEUK_EASYPREDICT reports it.
    output:
    path "${db_dir}", emit: database

    when:
    task.ext.when == null || task.ext.when

    script:
    // Some db names contain "/" (UniProtKB/Swiss-Prot), so db_dir is a filesystem-safe name,
    // distinct per db_name so storeDir caches never collide.
    db_dir = db_name.replaceAll('[^A-Za-z0-9_.-]', '_')
    def args = task.ext.args ?: ''
    """
    mkdir -p ${db_dir}

    metaeuk databases \\
        ${db_name} \\
        ${db_dir}/db \\
        tmp/ \\
        ${args} \\
        --threads ${task.cpus}

    rm -rf tmp/
    """

    stub:
    db_dir = db_name.replaceAll('[^A-Za-z0-9_.-]', '_')
    """
    mkdir -p ${db_dir}
    touch ${db_dir}/db.version
    touch ${db_dir}/db.dbtype
    """
}
