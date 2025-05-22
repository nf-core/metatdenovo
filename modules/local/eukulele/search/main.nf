process EUKULELE_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::eukulele=2.0.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:2.0.5--pyh723bec7_0' :
        'biocontainers/eukulele:2.0.5--pyh723bec7_0' }"

    input:
    tuple val(meta), path(fasta), val(dbname), path(eukdb)

    output:
    tuple val(meta), path("*/taxonomy_estimation/*.out.gz"), val("${dbname}") , emit: taxonomy_estimation
    tuple val(meta), path("*/taxonomy_counts/*.csv.gz")                       , emit: taxonomy_counts    , optional: true
    tuple val(meta), path("*/mets_full/diamond/*")                            , emit: diamond

    path "versions.yml"                                                       , emit: versions

    script:
    def args     = task.ext.args ?: ''
    def prefix   = task.ext.prefix ?: ("${dbname}" ? "${meta.id}_${dbname}" : "${meta.id}")
    def gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ./contigs/proteins.faa" : "mv ${fasta} contigs/proteins.faa"
    def database = dbname ? "--database ${dbname}" : ''

    """
    rc=0
    mkdir contigs
    $gunzip
    EUKulele \\
        $args \\
        $database \\
        --protein_extension .faa \\
        --reference_dir $eukdb \\
        -o ${prefix} \\
        --CPUs ${task.cpus} \\
        -s \\
        contigs || rc=\$?

    gzip ${prefix}/mets_full/diamond/*.out
    find ${prefix}/ -name "*.csv" | xargs gzip
    gzip ${prefix}/taxonomy_estimation/*.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukulele: \$(echo \$(EUKulele --version 2>&1) | sed -n 's/.* \\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS

    if [ \$rc -le 1 ]; then
        exit 0
    else
        exit \$rc;
    fi
    """
}
