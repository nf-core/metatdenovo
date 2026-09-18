process INTERPROSCAN_STRIPASTERISK {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(faa)

    output:
    tuple val(meta), path("*.faa"), emit: faa

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def is_compressed  = faa.getExtension() == "gz"
    def decompress     = is_compressed ? "gzip -c -d ${faa}" : "cat ${faa}"

    // ORF callers (e.g. Prodigal) mark stop codons with a trailing '*' -- InterProScan hard
    // rejects any '*' in a submitted sequence ("not a valid IUPAC amino acid character"), so
    // strip it here rather than in the vendored INTERPROSCAN module itself.
    """
    $decompress | sed '/^>/! s/\\*//g' > ${prefix}.faa
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.faa
    """
}
