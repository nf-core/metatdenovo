process KOFAMSCAN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kofamscan:1.3.0--hdfd78af_2':
        'quay.io/biocontainers/kofamscan:1.3.0--hdfd78af_2' }"

    input:
    tuple val(meta), path(fasta)
    path profiles
    path ko_list
    val  ko_list_url

    output:
    tuple val(meta), path('*.txt.gz'), optional: true, emit: txt
    tuple val(meta), path('*.tsv.gz'), optional: true, emit: tsv
    tuple val("${task.process}"), val('kofamscan'), eval("exec_annotation --version 2>&1 | sed 's/exec_annotation //;'"), topic: versions, emit: versions_kofamscan
    tuple val("${task.process}"), val('kofamscan_db'), val(kofamscan_db_provenance), topic: versions, emit: versions_kofamscan_db

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--format detail-tsv") ? "tsv" :
                    "txt"
    def is_compressed = fasta.toString().endsWith('.gz')
    def fasta_name     = is_compressed ? fasta.toString() - '.gz' : "$fasta"
    def gunzip         = is_compressed ? "gunzip -c ${fasta} > ${fasta_name}" : ''
    def cleanup        = is_compressed ? "rm -f ${fasta_name}" : ''
    // KEGG's FTP mirror has no release/date field of its own -- genome.jp always serves "latest"
    // at a stable path -- so record what IS knowable: the source URL and the date this local
    // copy was obtained, from ko_list's mtime. Read here (not in KOFAMSCAN_DOWNLOAD, and not with
    // eval()): ko_list is an already-completed upstream task's output, so it's guaranteed
    // materialized by the time this preamble runs, unlike reading it inside the same process
    // that's still staging a remote file(url) input, or a storeDir-cached process that may never
    // actually run its script again.
    kofamscan_db_provenance = "source=${ko_list_url}; retrieved=${new Date(ko_list.lastModified()).format('yyyy-MM-dd')}"
    """
    $gunzip

    exec_annotation \\
        -p $profiles \\
        -k $ko_list \\
        $args \\
        --cpu $task.cpus \\
        -o ${prefix}.${extension} \\
        ${fasta_name}

    gzip ${prefix}.${extension}

    $cleanup
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--format detail-tsv") ? "tsv" :
                    "txt"
    kofamscan_db_provenance = "source=${ko_list_url}; retrieved=${new Date(ko_list.lastModified()).format('yyyy-MM-dd')}"
    """
    touch ${prefix}.${extension}

    gzip ${prefix}.${extension}
    """
}
