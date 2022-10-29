process EUKULELE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::eukulele=2.0.2-0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:2.0.2--pyh723bec7_0' :
        'quay.io/biocontainers/eukulele:2.0.2--pyh723bec7_0' }"

    input:
    tuple val(meta), path(contigs_fasta)
    path(db)

    output:
    tuple val(meta), path("${meta.id}/taxonomy_estimation/*.out")                  , emit: taxonomy_estimation
    tuple val(meta), path("${meta.id}/taxonomy_counts/${meta.id}_all_*_counts.csv"), emit: taxonomy_counts
    tuple val(meta), path("${meta.id}/mets_full/diamond/*")                        , emit: diamond

    path "versions.yml"                                                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    rc=0
    mkdir contigs
    ln -s ${contigs_fasta} contigs

    EUKulele \\
        $args \\
        --reference_dir $db \\
        -o ${meta.id} \\
        --CPUs ${task.cpus} \\
        -s \\
        contigs || rc=\$?

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukulele: \$(echo \$(EUKulele --version 2>&1) | sed 's/Running EUKulele with command line arguments, as no valid configuration file was provided.//; s/The current EUKulele version is//g')
    END_VERSIONS

    if [ \$rc -le 1 ]; then
        exit 0
    else
        exit \$rc;
    fi
    """
}
