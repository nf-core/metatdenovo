process FORMAT_CLUSTERREPS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(clusters)

    output:
    tuple val(meta), path("${prefix}.representatives.txt"), emit: representatives
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    // MMseqs2's cluster table is one row per member, naming that member's representative in column
    // 1, so the distinct values of column 1 are the cluster representatives. Deduplicated in awk
    // rather than with sort -u so the list keeps the order the clusters appear in, and so this stays
    // a single streaming pass over a table that has one row per called gene on real data.
    """
    awk 'BEGIN { FS = "\\t" }
        !seen[\$1]++ { print \$1 }' ${clusters} \\
        > ${prefix}.representatives.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.representatives.txt
    """
}
