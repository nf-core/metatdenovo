process EGGNOG_DOWNLOAD {
    tag 'EggNOG'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    // Bundles awscli alongside eggnog-mapper: this process' storeDir can point at an s3:// path
    // (e.g. --eggnog_dbpath in conf/test_full.config), and Nextflow's AWS Batch executor stages
    // storeDir output to S3 by shelling out to `aws` from inside the task's own container --
    // which the plain eggnog-mapper image doesn't have, addresses #471.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/eggnog-mapper_awscli:34a6ca5baa89f396':
        'community.wave.seqera.io/library/eggnog-mapper_awscli:635add8f85922662' }"

    output:
    path "eggnog.db"                  , emit: eggnog_db
    path "eggnog_proteins.dmnd"       , emit: dmnd
    path "eggnog.taxa.db"             , emit: taxa_db
    path "eggnog.taxa.db.traverse.pkl", emit: pkl
    path "*"                          , emit: all

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # This commented for the moment since the tool tries to access a domain that doesn't exist anymore
    #download_eggnog_data.py \\
    #    $args \\
    #    -y \\
    #    --data_dir .

    # Temporary solution, until version 3 of the tool
    wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz
    gunzip eggnog.db.gz
    wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
    gunzip eggnog_proteins.dmnd.gz
    wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
    tar xzf eggnog.taxa.tar.gz
    """

    stub:
    """
    touch eggnog.db
    touch eggnog_proteins.dmnd
    touch eggnog.taxa.db
    touch eggnog.taxa.db.traverse.pkl
    """
}
