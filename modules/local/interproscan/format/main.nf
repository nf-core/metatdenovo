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

    // INTERPROSCAN is always run with --iprlookup --goterms (see conf/modules.config), which
    // fixes this column set; the raw TSV itself has no header row.
    """
    (
        printf 'orf\\tmd5\\tlength\\tanalysis\\tsignature_accession\\tsignature_description\\tstart\\tstop\\tscore\\tstatus\\tdate\\tinterpro_accession\\tinterpro_description\\tgo_terms\\n'
        cat ${tsv}
    ) | gzip -c > ${prefix}.interproscan.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    gzip -c /dev/null > ${prefix}.interproscan.tsv.gz
    """
}
