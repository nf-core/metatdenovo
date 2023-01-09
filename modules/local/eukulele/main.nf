process EUKULELE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::eukulele=2.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:2.0.3--pyh723bec7_0' :
        'quay.io/biocontainers/eukulele:2.0.3--pyh723bec7_0' }"

    input:
    tuple val(meta)  , path(fasta)
    tuple path(eukdb), val(namedb)

    output:
    tuple val(meta), path("${meta.id}_${namedb}/taxonomy_estimation/*.out")                            , emit: taxonomy_estimation
    tuple val(meta), path("${meta.id}_${namedb}/taxonomy_counts/${meta.id}_${namedb}_all_*_counts.csv"), emit: taxonomy_counts
    tuple val(meta), path("${meta.id}_${namedb}/mets_full/diamond/*")                                  , emit: diamond

    path "versions.yml"                                                                                , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input    = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input}" : ""

    """
    
    $gunzip

    rc=0
    mkdir contigs
    cp $input ./contigs/

    EUKulele \\
        $args \\
        --database $namedb \\
        --reference_dir $eukdb \\
        -o ${meta.id}_$namedb \\
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
