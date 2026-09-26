process KOFAMSCAN_DOWNLOAD {
    tag "KEGG data"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    // An s3:// storeDir is staged by aws inside this container, which the wget image lacks.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/bc/bceb5c307eb199ae3eca0eecfd71a0a4a918ce90e6fd96ecc749603426337823/data' :
        'community.wave.seqera.io/library/wget_awscli_gzip_tar:1fad694ee6322b7d' }"

    input:
    path ko_list_gz
    path profiles_targz

    output:
    path "ko_list"     , emit: ko_list
    path "profiles"    , emit: koprofiles

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    gunzip -c ${ko_list_gz} > ko_list

    tar -zxf ${profiles_targz}
    """

    stub:

    """
    touch ko_list
    mkdir profiles
    """
}
