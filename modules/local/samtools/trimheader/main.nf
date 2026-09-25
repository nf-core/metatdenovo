process SAMTOOLS_TRIMHEADER {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e9/e9b1a71ab018f22bcae5f800c93b0a6627516210c2f3bcf6c5a5af8cfdfb1d99/data'
        : 'community.wave.seqera.io/library/gawk_htslib_samtools:908360748e8474ae'}"

    input:
    tuple val(meta), path(bam, stageAs: 'input/*'), path(idxstats)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: bam.baseName
    """
    {
        samtools view -H ${bam} \\
            | awk -F '\\t' 'NR == FNR { if (\$3 + \$4 > 0) keep[\$1]; next } !/^@SQ/ || (substr(\$2, 4) in keep)' ${idxstats} -
        samtools view -@ ${task.cpus} ${bam}
    } | samtools view -b -@ ${task.cpus} -o ${prefix}.bam -
    """

    stub:
    def prefix = task.ext.prefix ?: bam.baseName
    """
    touch ${prefix}.bam
    """
}
