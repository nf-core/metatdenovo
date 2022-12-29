process KHMER_NORMALIZEBYMEDIAN {
    tag "${name}"
    label 'process_long'

    conda (params.enable_conda ? "bioconda::khmer=3.0.0a3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2' :
        'quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2' }"

    input:
    path pe_reads
    path se_reads
    val  name

    output:
    path "${name}.nm.fastq.gz", emit: reads
    path "${name}.nm.kh"      , emit: graph
    path "${name}.nm.log"     , emit: log
    path "versions.yml"       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ''
    pe_args = pe_reads ? "--paired" : ""
    se_args = se_reads ? "--unpaired-reads ${se_reads}" : ""
    files   = pe_reads ? pe_reads : se_reads

    """
    normalize-by-median.py \\
        -M ${task.memory.toGiga()}e9 \\
        --gzip $args \\
        --savegraph ${name}.nm.kh \\
        -o ${name}.nm.fastq.gz \\
        $args \\
        ${pe_args} \\
        ${se_args} \\
        ${files} 2>&1 | tee ${name}.nm.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        khmer_normalizebymedian: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
