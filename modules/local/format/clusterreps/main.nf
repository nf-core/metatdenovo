process FORMAT_CLUSTERREPS {
    tag "$meta.id"
    // Buffers every cluster's members, so memory scales with the number of loci.
    label 'process_medium'

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

    // MMseqs2 picks representatives by input order, so an unrelated upstream change can rename a
    // cluster. Re-pick the smallest member and sort, so ids depend on cluster content only.
    """
    awk 'BEGIN { FS = OFS = "\\t" }
        {
            if (!(\$1 in group_seen)) { group_seen[\$1] = 1; group_order[++n_groups] = \$1 }
            # if/else, not ternary: mawk creates the element before evaluating the RHS.
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
