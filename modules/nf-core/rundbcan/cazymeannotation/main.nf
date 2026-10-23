process RUNDBCAN_CAZYMEANNOTATION {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dbcan:5.2.9--pyhdfd78af_0' :
        'quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input_raw_data)
    path dbcan_db

    output:
    tuple val(meta), path("${prefix}_overview.tsv"), emit: cazyme_annotation
    tuple val(meta), path("${prefix}_dbCAN_hmm_results.tsv"), emit: dbcanhmm_results
    tuple val(meta), path("${prefix}_dbCANsub_hmm_results.tsv"), emit: dbcansub_results
    tuple val(meta), path("${prefix}_diamond.out"), emit: dbcandiamond_results
    tuple val("${task.process}"), val('rundbcan'), eval("run_dbcan version | sed 's/dbCAN version: //g'"), emit: versions_rundbcan, topic: versions
    tuple val("${task.process}"), val('dbcan_db'), val(dbcan_db_provenance), emit: versions_dbcan_db, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // RUNDBCAN_DATABASE's own output has no version/release tag (it's always "latest" from
    // dbCAN's built-in default mirror), so record what IS knowable: the date this local copy
    // was obtained, from the directory's mtime -- read here, not with eval(), since dbcan_db is
    // storeDir-cached and may never actually run its script again.
    dbcan_db_provenance = "source=run_dbcan database (built-in default mirror); retrieved=${new Date(dbcan_db.lastModified()).format('yyyy-MM-dd')}"
    """
    run_dbcan CAZyme_annotation \\
        --mode protein \\
        --db_dir ${dbcan_db} \\
        --input_raw_data ${input_raw_data} \\
        --output_dir . \\
        --threads ${task.cpus} \\
        ${args}

    mv overview.tsv ${prefix}_overview.tsv
    mv dbCAN_hmm_results.tsv ${prefix}_dbCAN_hmm_results.tsv
    mv dbCANsub_hmm_results.tsv ${prefix}_dbCANsub_hmm_results.tsv
    mv diamond.out ${prefix}_diamond.out
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    dbcan_db_provenance = "source=run_dbcan database (built-in default mirror); retrieved=${new Date(dbcan_db.lastModified()).format('yyyy-MM-dd')}"
    """
    touch ${prefix}_overview.tsv
    touch ${prefix}_dbCAN_hmm_results.tsv
    touch ${prefix}_dbCANsub_hmm_results.tsv
    touch ${prefix}_diamond.out
    """
}
