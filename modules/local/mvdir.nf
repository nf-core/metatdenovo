process MV_DIR {
    input:
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path('contigs'), emit: contigs_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir contigs
    mv ${fasta_file} contigs
    """
}
