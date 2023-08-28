process EUKULELE_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::eukulele=2.0.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:2.0.3--pyh723bec7_0' :
        'biocontainers/eukulele:2.0.3--pyh723bec7_0' }"

    input:
    tuple val(meta), path(fasta), val(dbname), path(eukdb)

    output:
    tuple val(meta), path("${meta.id}_*/taxonomy_estimation/*.out"), val("${db}")    , emit: taxonomy_estimation
    tuple val(meta), path("${meta.id}_*/taxonomy_counts/*_counts.csv")               , emit: taxonomy_counts
    tuple val(meta), path("${meta.id}_*/mets_full/diamond/*")                        , emit: diamond

    path "versions.yml"                                                              , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    input    = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ./contigs/${input}" : ""
    def database = dbname ? "--database ${dbname}" : ''
    db       = dbname ? "${dbname}" : 'default'

    """
    rc=0
    mkdir contigs
    $gunzip


    EUKulele \\
        $args \\
        $database \\
        --reference_dir $eukdb \\
        -o ${meta.id}_${db} \\
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
