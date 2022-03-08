process MV_DIR {
    input:
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path('contigs'), emit: contigs_dir

    script:
    """
    mkdir contigs
    mv ${fasta_file} contigs
    """
}
