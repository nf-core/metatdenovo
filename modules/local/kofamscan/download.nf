process KOFAMSCAN_DOWNLOAD {
    tag "KEGG data"
    label 'process_long'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hed695b0_4':
        'quay.io/biocontainers/gnu-wget:1.18--hed695b0_4' }"

    input:
    path(ko_list_file)
    path(koprofiles_dir)

    output:
    path("profiles/*"), emit: profiles
    path("ko_list/*") , emit: ko_list
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
    wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
    tar -xfz profiles.tar

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: 1.18
    END_VERSIONS
    """
}
