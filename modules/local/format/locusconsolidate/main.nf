process FORMAT_LOCUSCONSOLIDATE {
    tag "$meta.id"
    // Buffers provenance and members for every locus until END, so memory scales with the number
    // of called ORFs rather than being streamed; process_medium rather than process_low for that.
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(sorted_bed)

    output:
    tuple val(meta), path("${prefix}.gff.gz")           , emit: gff
    tuple val(meta), path("${prefix}.provenance.tsv.gz"), emit: provenance
    tuple val(meta), path("${prefix}.members.tsv.gz")   , emit: members
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"

    // Groups overlapping CDS calls into loci, from a BED sorted by (chrom, start) whose name column
    // is "<caller>:<original ID>" (see FORMAT_GFF2BED).
    //
    // The grouping rule is the point of this module: two calls join the same locus only if they
    // overlap on the same strand AND come from callers not already in that group. Plain interval
    // merging cannot express that second condition, and getting it wrong is not a corner case --
    // prokaryotic genes overlap each other routinely (the classic 4 bp ATGA operon overlap), so
    // merging on overlap alone fuses adjacent genes from ONE caller into a single locus. Measured on
    // this pipeline's own default test data, that fused 515 of 3470 Prodigal ORFs (14.8%) into 238
    // shared loci, silently summing the counts of genes that are not the same gene.
    //
    // With the caller condition, a single-caller run can never merge anything, so every locus has one
    // contributor, inherits that ORF's own ID, and the consolidated counts table is identical in
    // content to that caller's own -- the graceful-degradation property this level is supposed to
    // have. A locus with several contributors gets a fresh coordinate-derived ID instead.
    //
    // Grouping is a single greedy sweep in sorted order, so it is deterministic but not "optimal":
    // in a chain A(caller1) - B(caller2) - C(caller1), B binds to A and C starts a new locus.
    //
    // Because the sweep is greedy, WHICH call a locus absorbs depends on the order equal-start calls
    // arrive in, so the input is re-sorted here on a total order over the fields that matter --
    // chrom, start, end, strand, name -- rather than trusting the upstream sort. Ordering by
    // (chrom, start) alone leaves ties to be broken by whatever order the per-caller BEDs happened to
    // be concatenated in, which changes when a caller's task re-runs into a different work directory.
    // Measured on the full-size test data (91863 calls, 16335 coordinate groups holding more than one
    // call): 68829 loci with one caller's BED first, 68826 with the other's, 68828 shuffled. With the
    // sort below, all three give 68829. LC_ALL=C so the name comparison is byte order everywhere
    // rather than whatever collation the host locale supplies.
    //
    // Provenance and members are accumulated by locus ID and written at END rather than per group,
    // because a multi-exon gene contributes several non-overlapping groups that all inherit the same
    // ID. Emitting per group would repeat that ID with per-segment counts, and any consumer joining
    // on it would fan out and duplicate the locus's counts.
    """
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 ${sorted_bed} \\
        | awk 'BEGIN { FS = OFS = "\\t"; SEP = SUBSEP }

        # Close the open group for one contig/strand: assign the locus its ID and record provenance.
        function flush(key,   i, n, parts, cnt, id, sep, seen_member) {
            if (!(key in g_end)) return
            n = split(g_members[key], parts, SEP)
            delete seen_member
            cnt = 0
            for (i = 1; i <= n; i++) {
                if (!(parts[i] in seen_member)) { seen_member[parts[i]] = 1; cnt++ }
            }
            if (cnt == 1) {
                sep = index(parts[1], ":")
                id  = substr(parts[1], sep + 1)
            } else {
                id = "locus_" g_chrom[key] "_" (g_start[key] + 1) "_" g_end[key] "_" g_strand[key]
            }
            print g_chrom[key], "locus_consolidate", "CDS", g_start[key] + 1, g_end[key], ".", g_strand[key], ".", "ID=" id \\
                | "sort -k1,1 -k4,4n -k5,5n -k7,7 | gzip -c > ${prefix}.gff.gz"
            if (!(id in prov_seen)) { prov_seen[id] = 1; prov_order[++n_prov] = id }
            for (i = 1; i <= n; i++) {
                if (!((id SEP parts[i]) in member_seen)) {
                    member_seen[id SEP parts[i]] = 1
                    # if/else, not a ternary on the right of the assignment: some awks create the
                    # target array element before evaluating the right-hand side, which would make
                    # "id in prov_members" true on the first append and prepend an empty member.
                    if (id in prov_members) prov_members[id] = prov_members[id] SEP parts[i]
                    else                    prov_members[id] = parts[i]
                }
            }
            delete g_end[key]
        }

        {
            name   = \$4
            strand = \$6
            key    = \$1 SEP strand
            sep    = index(name, ":")
            caller = substr(name, 1, sep - 1)

            if ((key in g_end) && \$2 <= g_end[key] && index(g_callers[key], SEP caller SEP) == 0) {
                if (\$3 > g_end[key]) g_end[key] = \$3
                g_callers[key] = g_callers[key] caller SEP
                g_members[key] = g_members[key] SEP name
            } else {
                flush(key)
                g_chrom[key]   = \$1
                g_strand[key]  = strand
                g_start[key]   = \$2
                g_end[key]     = \$3
                g_callers[key] = SEP caller SEP
                g_members[key] = name
            }
        }

        END {
            for (key in g_end) flush(key)

            print "ID", "callers", "n_calls" | "gzip -c > ${prefix}.provenance.tsv.gz"
            print "ID", "caller", "orf"      | "gzip -c > ${prefix}.members.tsv.gz"
            for (p = 1; p <= n_prov; p++) {
                id = prov_order[p]
                n  = split(prov_members[id], parts, SEP)
                # Members and callers are insertion-sorted, and the GFF is piped through sort above,
                # so every output is a function of the locus content alone. Nothing here depends on
                # the order intervals happened to arrive in, which keeps the files reproducible even
                # if the upstream sort is not stable for ties.
                delete sorted_m
                for (i = 1; i <= n; i++) {
                    for (j = i - 1; j >= 1 && sorted_m[j] > parts[i]; j--) sorted_m[j + 1] = sorted_m[j]
                    sorted_m[j + 1] = parts[i]
                }
                n_callers = 0
                for (i = 1; i <= n; i++) {
                    sep    = index(sorted_m[i], ":")
                    caller = substr(sorted_m[i], 1, sep - 1)
                    orf    = substr(sorted_m[i], sep + 1)
                    print id, caller, orf | "gzip -c > ${prefix}.members.tsv.gz"
                    dup = 0
                    for (j = 1; j <= n_callers; j++) if (sorted[j] == caller) { dup = 1; break }
                    if (!dup) {
                        for (j = n_callers; j >= 1 && sorted[j] > caller; j--) sorted[j + 1] = sorted[j]
                        sorted[j + 1] = caller
                        n_callers++
                    }
                }
                callers = ""
                for (j = 1; j <= n_callers; j++) callers = (j == 1 ? sorted[j] : callers "," sorted[j])
                # n_calls counts the independent calls merged into this locus, i.e. distinct
                # contributing ORFs -- not exon segments, and not intervals.
                print id, callers, n | "gzip -c > ${prefix}.provenance.tsv.gz"
            }
        }'
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gzip -c /dev/null > ${prefix}.gff.gz
    gzip -c /dev/null > ${prefix}.provenance.tsv.gz
    gzip -c /dev/null > ${prefix}.members.tsv.gz
    """
}
