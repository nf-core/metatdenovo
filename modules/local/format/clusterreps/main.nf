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
    tuple val(meta), path("${prefix}.clusters.tsv")       , emit: clusters
    tuple val(meta), path("${prefix}.representatives.txt"), emit: representatives
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    // MMseqs2 names one member of each cluster as its representative, in column 1 of a table with one
    // row per member. Which member it picks depends on the order sequences were fed in, so an
    // unrelated change upstream -- a different sort, an extra caller -- can silently rename a cluster
    // even though its membership is identical. That name is the cluster id in the counts table and
    // the sequence id handed to the annotation tools, so it should be a function of the cluster
    // content and nothing else.
    //
    // Re-pick it as the lexicographically smallest member and rewrite the table, so the counts table
    // and the annotated representatives can never disagree about what a cluster is called. Clusters
    // are emitted in order of first appearance, and members sorted within a cluster, so the output
    // does not depend on input order either.
    """
    awk 'BEGIN { FS = OFS = "\\t" }
        {
            if (!(\$1 in group_seen)) { group_seen[\$1] = 1; group_order[++n_groups] = \$1 }
            # if/else rather than a ternary on the right of the assignment: some awks create the
            # target element before evaluating the right-hand side, which would prepend an empty
            # member on the first append.
            if (\$1 in members) members[\$1] = members[\$1] SUBSEP \$2
            else                members[\$1] = \$2
        }
        END {
            for (g = 1; g <= n_groups; g++) {
                n = split(members[group_order[g]], parts, SUBSEP)
                delete sorted
                for (i = 1; i <= n; i++) {
                    for (j = i - 1; j >= 1 && sorted[j] > parts[i]; j--) sorted[j + 1] = sorted[j]
                    sorted[j + 1] = parts[i]
                }
                print sorted[1] > "${prefix}.representatives.txt"
                for (i = 1; i <= n; i++) print sorted[1], sorted[i] > "${prefix}.clusters.tsv"
            }
        }' ${clusters}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.clusters.tsv
    touch ${prefix}.representatives.txt
    """
}
