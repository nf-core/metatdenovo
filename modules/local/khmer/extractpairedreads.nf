process KHMER_EXTRACTPAIREDREADS {
    tag "$name"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::khmer=3.0.0a3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    path reads
    val  name

    output:
    path "${name}.ep.pe.fastq.gz", emit: pairs
    path "${name}.ep.se.fastq.gz", emit: singles
    path "${name}.ep.log"        , emit: log
    path "versions.yml"          , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ''
    
    """
    extract-paired-reads.py \\
        --gzip $args \\
        -p ${name}.ep.pe.fastq.gz \\
        -s ${name}.ep.se.fastq.gz \\
        ${reads} 2>&1 | tee ${name}.ep.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer_extractreadpairs: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
