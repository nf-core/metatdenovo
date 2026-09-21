process EGGNOG_DOWNLOAD {
    tag 'EggNOG'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // storeDir can point at an s3:// path; Nextflow's AWS Batch executor stages that copy by
    // shelling out to `aws` inside the task's own container, which the plain eggnog-mapper image lacks.
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/eggnog-mapper_awscli:34a6ca5baa89f396':
        'community.wave.seqera.io/library/eggnog-mapper_awscli:635add8f85922662' }"

    input:
    path eggnog_db_gz
    path eggnog_dmnd_gz
    path eggnog_taxa_targz

    output:
    path "eggnog.db"                  , emit: eggnog_db
    path "eggnog_proteins.dmnd"       , emit: dmnd
    path "eggnog.taxa.db"             , emit: taxa_db
    path "eggnog.taxa.db.traverse.pkl", emit: pkl
    path "*"                          , emit: all

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Temporary solution, until version 3 of the tool: download_eggnog_data.py
    # itself tries to reach a domain that no longer exists.
    gunzip -c ${eggnog_db_gz} > eggnog.db
    gunzip -c ${eggnog_dmnd_gz} > eggnog_proteins.dmnd
    tar xzf ${eggnog_taxa_targz}
    """

    stub:
    """
    touch eggnog.db
    touch eggnog_proteins.dmnd
    touch eggnog.taxa.db
    touch eggnog.taxa.db.traverse.pkl
    """
}
