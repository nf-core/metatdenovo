process KOFAMSCAN_DOWNLOAD {
    tag "KEGG data"
    label 'process_long'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--h36e9172_9':
        'quay.io/biocontainers/gnu-wget:1.18--h36e9172_9' }"

    output:
    path "ko_list"     , emit: ko_list
    path "profiles"    , emit: koprofiles
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
    gunzip ko_list.gz

    wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
    tar -zxf profiles.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | grep 'GNU Wget' | sed 's/GNU Wget \\([0-9.]\\+\\) .*/\\1/')
    END_VERSIONS
    """
}
