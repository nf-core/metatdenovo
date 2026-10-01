process FORMAT_LOCUSFAA {
    tag "$meta.id"
    // Holds all protein sequences in memory, so memory scales with total protein length.
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(members), val(callers), path(faas)

    output:
    tuple val(meta), path("${prefix}.locus.faa.gz"), emit: faa
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    // Ids are rewritten to "<caller>:<ORF id>", the members-table key; only the first header token is the id.
    def normalise = [callers, faas].transpose().collect { caller, faa ->
        def cat_input = faa.name.endsWith('.gz') ? "gunzip -c ${faa}" : "cat ${faa}"
        "${cat_input} | awk -v caller='${caller}' '/^>/ { split(substr(\$0, 2), a, \" \"); print \">\" caller \":\" a[1]; next } { print }' >> all_proteins.faa"
    }.join('\n    ')

    // Longest sequence per locus, ties broken on smallest "<caller>:<ORF id>";
    // emitted in members-table order since awk for-in order is unspecified.
    """
    : > all_proteins.faa
    ${normalise}

    gunzip -c ${members} \\
        | awk 'BEGIN { FS = OFS = "\\t" }
            NR == FNR {
                if (FNR > 1) {
                    locus[\$2 ":" \$3] = \$1
                    if (!(\$1 in seen_locus)) {
                        seen_locus[\$1] = 1
                        locus_order[++n_loci] = \$1
                    }
                }
                next
            }
            /^>/ {
                key = substr(\$0, 2)
                key_order[++n_keys] = key
                next
            }
            { seq[key] = seq[key] \$0 }
            END {
                for (i = 1; i <= n_keys; i++) {
                    k = key_order[i]
                    if (!(k in locus)) continue
                    l = locus[k]
                    s = seq[k]
                    if (!(l in best_key) \\
                        || length(s) > length(best_seq[l]) \\
                        || (length(s) == length(best_seq[l]) && k < best_key[l])) {
                        best_key[l] = k
                        best_seq[l] = s
                    }
                }
                n_emitted = 0
                for (j = 1; j <= n_loci; j++) {
                    l = locus_order[j]
                    if (l in best_seq) { print ">" l "\\n" best_seq[l]; n_emitted++ }
                    else if (n_missing++ < 5) missing = missing " " l
                }
                # A locus without a protein would silently lose its counts downstream, so fail.
                if (n_emitted < n_loci) {
                    printf "ERROR: %d of %d loci had no matching protein sequence, e.g.%s\\n", \\
                        n_loci - n_emitted, n_loci, missing > "/dev/stderr"
                    exit 1
                }
            }' - all_proteins.faa \\
        | gzip -c > ${prefix}.locus.faa.gz

    rm all_proteins.faa
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -c /dev/null > ${prefix}.locus.faa.gz
    """
}
