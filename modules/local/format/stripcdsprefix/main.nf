process FORMAT_STRIPCDSPREFIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11':
        'biocontainers/gzip:1.11' }"

    input:
    tuple val(meta), path(counts)

    output:
    tuple val(meta), path("${prefix}.counts.tsv.gz"), emit: counts
    tuple val("${task.process}"), val('gzip'), eval('gzip --version  2>&1 | grep "^gzip" | sed "s/^gzip \\([0-9.]\\+\\).*/\\1/"'), emit: versions_gzip, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # TransDecoder prefixes gff ORF IDs with "cds.", but not fasta ones; orf is the first column.
    # Output name equals the staged input symlink: write elsewhere and rename.
    gzip -dc ${counts} | sed 's/^cds\\.//' | gzip -c > stripped.counts.tsv.gz
    mv -f stripped.counts.tsv.gz ${prefix}.counts.tsv.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > stub.counts.tsv.gz
    mv -f stub.counts.tsv.gz ${prefix}.counts.tsv.gz
    """
}
