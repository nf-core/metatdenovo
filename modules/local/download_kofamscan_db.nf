process DOWNLOAD_KODB {
    tag "downlaod_KO_databases"
    label 'process_medium'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hed695b0_4':
        'quay.io/biocontainers/gnu-wget:1.18--hed695b0_4' }"

    input:

    output:
    tuple val(meta), path("profiles/*"), emit: profiles
    tuple val(meta), path("ko_list/*") , emit: ko_list
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    wget https://www.genome.jp/ftp/db/kofam/ko_list.gz 
    wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
    gunzip ko_list.gz
    gunzip profiles.tar.gz
    tar -xf profiles.tar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.18
    END_VERSIONS
    """
}
