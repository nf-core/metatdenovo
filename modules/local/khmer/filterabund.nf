process KHMER_FILTERABUND {
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::khmer=3.0.0a3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    path graph
    path reads
    val  name

    output:
    path "${name}.fa.fastq.gz", emit: reads
    path "${name}.fa.log"     , emit: log
    path "versions.yml"       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ''

    """
    filter-abund.py \\
        --gzip \\
        -T ${task.cpus}   \\
        $args \\
        -o ${name}.fa.fastq.gz \\
        ${reads} \\
        ${graph} 2>&1 | tee ${name}.fa.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer_filterabund: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
