process FORMAT_LOCUS_FAA {
    tag "$meta.id"
    label 'process_low'

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

    // Each caller's protein fasta gets its sequence ids rewritten to the same "<caller>:<ORF id>"
    // form FORMAT_GFF2BED puts in the bed name column, which is what the members table is keyed on.
    // Only the first whitespace-delimited token of a header is the id -- Prokka appends a product
    // description and Prodigal a "# start # end # ..." coordinate block. Compression differs per
    // caller (Prokka's faa is gzipped, MetaEuk's fas is not), hence the per-file decompression
    // choice made here rather than in the shell.
    def normalise = [callers, faas].transpose().collect { caller, faa ->
        def cat_input = faa.name.endsWith('.gz') ? "gunzip -c ${faa}" : "cat ${faa}"
        "${cat_input} | awk -v caller='${caller}' '/^>/ { split(substr(\$0, 2), a, \" \"); print \">\" caller \":\" a[1]; next } { print }' >> all_proteins.faa"
    }.join('\n    ')

    // One protein per consolidated locus, named after the locus. A locus with a single contributor
    // trivially takes that contributor's sequence; a locus several callers agreed on takes the
    // longest of their near-identical sequences, breaking ties on the lexicographically smallest
    // "<caller>:<ORF id>" so the choice is reproducible rather than dependent on file order.
    // Emission follows the order loci first appear in the members table -- awk's for-in iteration
    // order is unspecified and would make the output non-deterministic between runs.
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
                for (j = 1; j <= n_loci; j++) {
                    l = locus_order[j]
                    if (l in best_seq) print ">" l "\\n" best_seq[l]
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
