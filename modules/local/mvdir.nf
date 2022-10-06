process MV_DIR {
    label 'process_single'

    input:
    tuple val(meta), path(fasta_file)

    output:
    tuple val(meta), path('contigs'), emit: contigs_dir
    file "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir contigs
    mv ${fasta_file} contigs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mv: no version
    END_VERSIONS
    """
}
