process KOFAMSCAN_DOWNLOAD {
    tag "KEGG data"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    // storeDir can point at an s3:// path; Nextflow's AWS Batch executor stages that copy by
    // shelling out to `aws` inside the task's own container, which the plain wget image lacks.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget_awscli:340260e7e9dd32f7':
        'community.wave.seqera.io/library/wget_awscli:9510e6a6af2abe94' }"

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
