process TRANSDECODER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/transdecoder:5.7.1--pl5321hdfd78af_0' :
    'biocontainers/transdecoder:5.7.1--pl5321hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}/*.pep") , emit: pep
    tuple val(meta), path("${prefix}/*.gff3"), emit: gff
    tuple val(meta), path("${prefix}/*.cds") , emit: cds
    tuple val(meta), path("${prefix}/*.bed") , emit: bed
    tuple val("${task.process}"), val('transdecoder'), eval('TransDecoder.LongOrfs --version | sed -e "s/TransDecoder.LongOrfs //g"'), emit: versions_transdecoder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    TransDecoder.LongOrfs \\
        $args \\
        -O $prefix \\
        -t \\
        $fasta

    TransDecoder.Predict \\
        $args \\
        -O $prefix \\
        -t \\
        $fasta
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}/empty.pep
    touch ${prefix}/empty.gff3
    touch ${prefix}/empty.cds
    touch ${prefix}/empty.bed
    """
}
