process SPADES {
    tag "$assembly"
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::spades=3.15.4' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/spades:3.15.4--h95f258a_0' :
        'quay.io/biocontainers/spades:3.15.4--h95f258a_0' }"

    input:
    tuple val(meta), path(illumina)
    val assembly

    output:
    tuple val(meta), path('*.scaffolds.fa.gz')    , optional:true, emit: scaffolds
    tuple val(meta), path('*.contigs.fa.gz')      , optional:true, emit: contigs
    path('*.transcripts.fa.gz')                   , optional:true, emit: transcripts
    tuple val(meta), path('*.gene_clusters.fa.gz'), optional:true, emit: gene_clusters
    tuple val(meta), path('*.assembly.gfa.gz')    , optional:true, emit: gfa
    tuple val(meta), path('*.log')                , emit: log
    path  "versions.yml"                          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def maxmem = task.memory.toGiga()

    """
    spades.py \\
        $args \\
        --threads $task.cpus \\
        --memory $maxmem \\
        --pe-12 1  --pe-12 2  \\
        --12 $illumina \\
        -o ./
    
    mv spades.log ${assembly}.spades.log

    if [ -f scaffolds.fasta ]; then
        mv scaffolds.fasta ${assembly}.scaffolds.fa
        gzip -n ${assembly}.scaffolds.fa
    fi
    if [ -f contigs.fasta ]; then
        mv contigs.fasta ${assembly}.contigs.fa
        gzip -n ${assembly}.contigs.fa
    fi
    if [ -f transcripts.fasta ]; then
        mv transcripts.fasta ${assembly}.transcripts.fa
        gzip -n ${assembly}.transcripts.fa
    fi
    if [ -f assembly_graph_with_scaffolds.gfa ]; then
        mv assembly_graph_with_scaffolds.gfa ${assembly}.assembly.gfa
        gzip -n ${assembly}.assembly.gfa
    fi
    if [ -f gene_clusters.fasta ]; then
        mv gene_clusters.fasta ${assembly}.gene_clusters.fa
        gzip -n ${assembly}.gene_clusters.fa
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed 's/^.*SPAdes genome assembler v//; s/ .*\$//')
    END_VERSIONS
    """
}
