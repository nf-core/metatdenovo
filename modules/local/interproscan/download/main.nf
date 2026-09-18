process INTERPROSCAN_DATABASE {
    tag "InterProScan data"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    // storeDir can point at an s3:// path; Nextflow's AWS Batch executor stages that copy by
    // shelling out to `aws` inside the task's own container, which the plain wget image lacks.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/wget_awscli:340260e7e9dd32f7':
        'community.wave.seqera.io/library/wget_awscli:9510e6a6af2abe94' }"

    input:
    val database_url

    output:
    path "interproscan_db", emit: db

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget ${database_url} -O interproscan_db.tar.gz

    # The upstream tarball extracts to interproscan-<version>/{bin,data,lib,...} -- the module
    # supplies its own bin/lib via the container's conda install, so only the (large) data/
    # subdirectory is needed, staged directly as INTERPROSCAN's `data` input.
    tar -zxf interproscan_db.tar.gz
    mv interproscan-*/data interproscan_db
    rm -rf interproscan-* interproscan_db.tar.gz
    """

    stub:
    """
    mkdir interproscan_db
    """
}
