process EGGNOG_MAPPER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.9--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path(eggnog_files), stageAs: 'eggnog/*'

    output:
    tuple val(meta), path("*.emapper.hits.gz")                , emit: hits
    tuple val(meta), path("*.emapper.seed_orthologs.gz")      , emit: seed_orthologs
    tuple val(meta), path("*.emapper.annotations.gz")         , emit: annotations
    tuple val(meta), path("*.emapper.tsv.gz")                 , emit: emappertsv
    tuple val(meta), path("*.emapper.annotations.xlsx")       , emit: xlsx,      optional: true
    tuple val(meta), path("*.emapper.orthologs.gz")           , emit: orthologs, optional: true
    tuple val(meta), path("*.emapper.genepred.fasta.gz")      , emit: genepred,  optional: true
    tuple val(meta), path("*.emapper.gff.gz")                 , emit: gff,       optional: true
    tuple val(meta), path("*.emapper.no_annotations.fasta.gz"), emit: no_anno,   optional: true
    tuple val(meta), path("*.emapper.pfam.gz")                , emit: pfam,      optional: true
    path "versions.yml"                                       , emit: versions

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    input    = fasta.name.endsWith(".gz") ? fasta.baseName : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input}" : ""

    """
    $gunzip

    emapper.py \\
        $args \\
        --cpu $task.cpus \\
        --data_dir eggnog \\
        --output $prefix \\
        -i $input

    gzip ${prefix}.emapper.*
    zgrep -v '^##' ${prefix}.emapper.annotations | \\
        sed 's/^#// ; /^query/s/.*/\\L&/ ; s/query/orf/' | \\
        gzip -c > ${prefix}.emapper.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//' | sed 's/ \\/ Expected.*//')
    END_VERSIONS
    """

    stub:
    """
    touch test.emapper.hits
    touch test.emapper.seed_orthologs
    touch test.emapper.annotations

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//' | sed 's/ \\/ Expected.*//')
    END_VERSIONS
    """
}
