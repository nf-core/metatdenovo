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
    path "${db_name}", emit: database

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    metaeuk databases \\
        ${db_name} \\
        ${db_name}/${db_name} \\
        tmp/ \\
        ${args} \\
        --threads ${task.cpus}

    rm -rf tmp/
    """

    stub:
    """
    mkdir -p ${db_name}
    touch ${db_name}/${db_name}.version
    touch ${db_name}/${db_name}.dbtype
    """
}
