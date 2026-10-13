process INTERPROSCAN_FORMAT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*.interproscan.tsv.gz"), emit: tsv
    tuple val("${task.process}"), val('gzip'), eval('gzip --version 2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"

    // InterProScan's raw TSV has no header row and is ragged: a match with no InterPro entry
    // omits the interpro_accession/interpro_description/go_terms/pathways columns entirely
    // rather than leaving them blank, so rows can be 11, 13, 14 or 15 columns wide depending on
    // what each match actually has, independent of which --iprlookup/--goterms/--pathways flags
    // were passed. Pad every row out to the full 15-column superset with '-' (InterProScan's own
    // missing-value marker) so the table is rectangular regardless of ext.args.
    """
    (
        printf 'orf\\tmd5\\tlength\\tanalysis\\tsignature_accession\\tsignature_description\\tstart\\tstop\\tscore\\tstatus\\tdate\\tinterpro_accession\\tinterpro_description\\tgo_terms\\tpathways\\n'
        awk 'BEGIN { FS = OFS = "\\t" } { for (i = NF + 1; i <= 15; i++) \$i = "-"; print }' ${tsv}
    ) | gzip -c > ${prefix}.interproscan.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gzip -c /dev/null > ${prefix}.interproscan.tsv.gz
    """
}
