process KOFAMSCAN_DOWNLOAD {
    tag "KEGG data"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3b/3b54fa9135194c72a18d00db6b399c03248103f87e43ca75e4b50d61179994b3/data':
        'community.wave.seqera.io/library/wget:1.21.4--8b0fcde81c17be5e' }"

    input:
    val ko_list_url
    val profiles_url

    output:
    path "ko_list"     , emit: ko_list
    path "profiles"    , emit: koprofiles

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    wget ${ko_list_url} -O ko_list.gz
    gunzip ko_list.gz

    wget ${profiles_url} -O profiles.tar.gz
    tar -zxf profiles.tar.gz
    """

    stub:

    """
    touch ko_list
    mkdir profiles
    """
}
