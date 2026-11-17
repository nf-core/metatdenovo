process EUKULELE_SEARCH {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:2.1.2--pyhdfd78af_0' :
        'biocontainers/eukulele:2.1.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), val(dbname), path(eukdb)

    output:
    tuple val(meta), path("${prefix}/taxonomy_estimation/*.out.gz"), val("${dbname}"), emit: taxonomy_estimation
    tuple val(meta), path("${prefix}/taxonomy_counts/*.csv.gz")                      , emit: taxonomy_counts, optional: true
    tuple val(meta), path("${prefix}/mets_full/diamond/*")                           , emit: diamond
    path("versions.yml"), emit: versions, topic: versions

    script:
    def args     = task.ext.args ?: ''
    prefix       = task.ext.prefix ?: ("${dbname}" ? "${meta.id}_${dbname}" : "${meta.id}")
    def gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ./contigs/proteins.faa" : "mv ${fasta} contigs/proteins.faa"
    def database = dbname ? "--database ${dbname}" : ''

    // BUSCO is not in the container, so without --no_busco every run exits 1.
    // EUKulele writes its errors only to ${prefix}/log/*.err, hence the tail on failure.
    """
    mkdir contigs
    $gunzip
    EUKulele \\
        $args \\
        $database \\
        --protein_extension .faa \\
        --reference_dir $eukdb \\
        -o ${prefix} \\
        --CPUs ${task.cpus} \\
        --no_busco \\
        -s \\
        contigs || { rc=\$?; tail -n 20 ${prefix}/log/*.err >&2 || true; exit \$rc; }

    gzip ${prefix}/*_full/diamond/*.out
    find ${prefix}/ -name "*.csv" -exec gzip {} +
    gzip ${prefix}/taxonomy_estimation/*.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukulele: \$(EUKulele --version 2>&1 | grep "current EUKulele version" | sed "s/The current EUKulele version is //")
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: ("${dbname}" ? "${meta.id}_${dbname}" : "${meta.id}")
    """
    mkdir -p ${prefix}/taxonomy_estimation ${prefix}/taxonomy_counts ${prefix}/mets_full/diamond
    gzip -c /dev/null > ${prefix}/taxonomy_estimation/empty.out.gz
    gzip -c /dev/null > ${prefix}/taxonomy_counts/empty.csv.gz
    gzip -c /dev/null > ${prefix}/mets_full/diamond/empty

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukulele: 2.1.2
    END_VERSIONS
    """
}
