process KOFAMSCAN_DOWNLOAD {
    tag "KEGG data"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data':
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"

    output:
    path "ko_list"     , emit: ko_list
    path "profiles"    , emit: koprofiles

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
    gunzip ko_list.gz

    wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
    tar -zxf profiles.tar.gz
    """
}
