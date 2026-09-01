process METAEUK_DOWNLOAD {
    tag "$db_name"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaeuk:6.a5d39d9--pl5321hf1761c0_2':
        'quay.io/biocontainers/metaeuk:6.a5d39d9--pl5321hf1761c0_2' }"

    input:
    val db_name

    // storeDir (below, in conf/modules.config) only allows val/path outputs, so this process
    // reports no version itself -- METAEUK_EASYPREDICT already reports the same tool's version
    // on every run that reaches it.
    output:
    path "${db_dir}", emit: database

    when:
    task.ext.when == null || task.ext.when

    script:
    // Some names `metaeuk databases -h` lists contain a `/` (e.g. UniProtKB/Swiss-Prot), which
    // can't be a directory/prefix name -- db_dir is a filesystem-safe stand-in used only for
    // that, kept distinct per db_name so storeDir doesn't reuse one database's cache for
    // another. The prefix inside db_dir doesn't need to match db_name: METAEUK_EASYPREDICT finds
    // it later from whichever *.version file is actually there.
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
