// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process KHMER_FILTERABUND {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::khmer=3.0.0a3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/khmer:3.0.0a3--py37haa7609a_2"
    } else {
        container "quay.io/biocontainers/khmer:3.0.0a3--py37haa7609a_2"
    }

    input:
    path graph
    path reads
    val  name

    output:
    path "${name}.fa.fastq.gz", emit: reads
    path "versions.yml"       , emit: versions

    script:
    
    """
    filter-abund.py \\
        --gzip \\
        -T ${task.cpus}   \\
        ${options.args} \\
        -o ${name}.fa.fastq.gz \\
        ${reads} \\
        ${graph}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
