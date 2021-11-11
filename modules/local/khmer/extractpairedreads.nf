// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from '../functions'

params.options = [:]
options        = initOptions(params.options)

process KHMER_EXTRACTPAIREDREADS {
    tag "$name"
    label 'process_long'
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
    path reads
    val  name

    output:
    path "${name}.ep.pe.fastq.gz", emit: pairs
    path "${name}.ep.se.fastq.gz", emit: singles
    path "${name}.ep.log"        , emit: log
    path "versions.yml"          , emit: versions

    script:
    
    """
    extract-paired-reads.py \\
        --gzip ${options.args} \\
        -p ${name}.ep.pe.fastq.gz \\
        -s ${name}.ep.se.fastq.gz \\
        ${reads} 2>&1 | tee ${name}.ep.log

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( normalize-by-median.py --version 2>&1 | grep ^khmer | sed 's/^khmer //' )
    END_VERSIONS
    """
}
