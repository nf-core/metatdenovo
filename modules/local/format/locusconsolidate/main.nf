process FORMAT_LOCUSCONSOLIDATE {
    tag "$meta.id"
    // Buffers every locus until END, so memory scales with ORF count.
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

    // Groups overlapping CDS calls (BED name "<caller>:<ID>") into loci. Calls join a locus only if
    // they overlap on the same strand and come from a caller not yet in it: same-caller genes overlap
    // routinely and must stay separate. Single-contributor loci keep the ORF ID.
    // The sweep is greedy, so the input is re-sorted on a total order to make ties deterministic.
    """
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k6,6 -k4,4 ${sorted_bed} \\
        | awk 'BEGIN { FS = OFS = "\\t"; SEP = SUBSEP }

        function flush(key,   i, n, parts, cnt, id, sep, seen_member, n_found) {
            if (!(key in g_end)) return
            n = split(g_members[key], parts, SEP)
            delete seen_member
            cnt = 0
            for (i = 1; i <= n; i++) {
                if (!(parts[i] in seen_member)) { seen_member[parts[i]] = 1; cnt++ }
            }
            # Reuse the locus a member ORF (an earlier exon) is already in, so a gene has one locus.
            id = ""
            n_found = 0
            for (i = 1; i <= n; i++) {
                if (parts[i] in orf_locus && orf_locus[parts[i]] != id) {
                    id = orf_locus[parts[i]]
                    n_found++
                }
            }
            if (n_found > 1) {
                printf "ERROR: group at %s:%d-%d(%s) holds ORFs from %d different loci, which cannot be merged here\\n", \\
                    g_chrom[key], g_start[key] + 1, g_end[key], g_strand[key], n_found > "/dev/stderr"
                exit 1
            }
            if (id == "") {
                if (cnt == 1) {
                    sep = index(parts[1], ":")
                    id  = substr(parts[1], sep + 1)
                } else {
                    id = "locus_" g_chrom[key] "_" (g_start[key] + 1) "_" g_end[key] "_" g_strand[key]
                }
            }
            for (i = 1; i <= n; i++) orf_locus[parts[i]] = id
            print g_chrom[key], "locus_consolidate", "CDS", g_start[key] + 1, g_end[key], ".", g_strand[key], ".", "ID=" id \\
                | "sort -k1,1 -k4,4n -k5,5n -k7,7 | gzip -c > ${prefix}.gff.gz"
            if (!(id in prov_seen)) { prov_seen[id] = 1; prov_order[++n_prov] = id }
            for (i = 1; i <= n; i++) {
                if (!((id SEP parts[i]) in member_seen)) {
                    member_seen[id SEP parts[i]] = 1
                    # if/else, not ternary: mawk creates the element before evaluating the RHS.
                    if (id in prov_members) prov_members[id] = prov_members[id] SEP parts[i]
                    else                    prov_members[id] = parts[i]
                }
            }
            # Drop all group state; leftovers accumulate one row per contig.
            delete g_chrom[key]
            delete g_start[key]
            delete g_end[key]
            delete g_strand[key]
            delete g_callers[key]
            delete g_members[key]
        }

        # Groups never span contigs; flushing per contig bounds memory and fixes emission order.
        function flush_contig(   i, n, keys) {
            n = 0
            for (i in open_keys) keys[++n] = i
            for (i = 1; i <= n; i++) flush(keys[i])
            delete open_keys
            delete orf_locus
            delete member_seen
        }

        {
            name   = \$4
            strand = \$6
            key    = \$1 SEP strand
            sep    = index(name, ":")
            caller = substr(name, 1, sep - 1)

            if (\$1 != contig) {
                if (contig != "") flush_contig()
                contig = \$1
            }
            open_keys[key] = 1

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
            if (contig != "") flush_contig()

            print "ID", "callers", "n_calls" | "gzip -c > ${prefix}.provenance.tsv.gz"
            print "ID", "caller", "orf"      | "gzip -c > ${prefix}.members.tsv.gz"
            for (p = 1; p <= n_prov; p++) {
                id = prov_order[p]
                n  = split(prov_members[id], parts, SEP)
                # Sort members and callers so output depends only on locus content.
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
                # n_calls: distinct contributing ORFs, not exon segments.
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
