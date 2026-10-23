process METAEUK_EASYPREDICT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/metaeuk:6.a5d39d9--pl5321hf1761c0_2':
        'quay.io/biocontainers/metaeuk:6.a5d39d9--pl5321hf1761c0_2' }"

    input:
    tuple val(meta), path(fasta)
    path(database)

    output:
    tuple val(meta), path("${prefix}.fas")      , emit: faa
    tuple val(meta), path("${prefix}.codon.fas"), emit: codon
    tuple val(meta), path("*.tsv")              , emit: tsv
    tuple val(meta), path("*.gff")              , emit: gff
    tuple val("${task.process}"), val('metaeuk'), eval("metaeuk | sed -n 's/.*Version: //p'"), emit: versions_metaeuk, topic: versions
    tuple val("${task.process}"), val('metaeuk_db'), val(db_version), emit: versions_metaeuk_db, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    // `database` is already staged on disk by the time this preamble runs, so its *.version file
    // (written by `metaeuk databases`, same file checked by workflows/metatdenovo.nf) can be read
    // here directly -- unlike `eval()`, this doesn't need the task to actually execute.
    def db_version_file = database.isDirectory() ? database.listFiles().find { f -> f.name.endsWith('.version') } : null
    db_version = db_version_file ? db_version_file.text.trim() : 'unknown'
    """
    if [ -d ${database} ]; then
        ## if supplying an mmseqs database as a directory, metaeuk requires the basename of the database
        DBBASE=`find ${database}/ -name "*.version" -exec sh -c 'file=\$(basename {}); echo \${file%%.*}' \\;`
        DB=`echo "${database}/\${DBBASE}"`
    else
        DB=${database}
    fi

    metaeuk easy-predict \\
        ${fasta} \\
        \${DB} \\
        ${prefix} \\
        tmp/ \\
        ${args} \\
        --threads ${task.cpus}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def db_version_file = database.isDirectory() ? database.listFiles().find { f -> f.name.endsWith('.version') } : null
    db_version = db_version_file ? db_version_file.text.trim() : 'unknown'
    """
    touch ${prefix}.fas
    touch ${prefix}.codon.fas
    touch ${prefix}.headersMap.tsv
    touch ${prefix}.gff
    """
}
