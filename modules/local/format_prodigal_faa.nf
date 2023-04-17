process FORMAT_PRODIGAL_FAA {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gzip=1.12" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzip:1.11 ':
        'quay.io/biocontainers/gzip:1.11 ' }"

    input:
    tuple val(meta), path (fasta)

    output:
    tuple val(meta), path("${prefix}_format.faa.gz"), emit: format_faa
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    input    = fasta =~ /\.gz$/ ? fasta.name.take(fasta.name.lastIndexOf('.')) : fasta
    gunzip   = fasta =~ /\.gz$/ ? "gunzip -c ${fasta} > ${input}" : ""

    """
    $gunzip
    sed 's/\\(\\(k[0-9]\\+_[0-9]\\+\\).*\\)ID=[0-9]\\+\\(_[0-9]\\+\\)/\\1ID=\\2\\3/g' $input > ${prefix}_format.faa
    sed -i 's/\\(\\(NODE_[0-9]\\+\\).*\\)ID=[0-9]\\+\\(_[0-9]\\+\\)/\\1ID=\\2\\3/g' ${prefix}_format.faa
    gzip ${prefix}_format.faa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gzip: \$( echo \$(gzip --version 2>&1)| sed 's/.* gzip//')
    END_VERSIONS
    """
}
