process EGGNOG_DOWNLOAD {
    tag 'EggNOG'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.9--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0' }"

    output:
    path "eggnog_data"           , emit: eggnog_data_dir
    path "eggnog_proteins.dmnd"  , emit: dmnd
    path "*"                     , emit: all

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
    mkdir eggnog_data
    wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.db.gz -O eggnog_data/eggnog.db.gz
    gunzip eggnog_data/eggnog.db.gz
    wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz
    gunzip eggnog_proteins.dmnd.gz
    wget http://eggnog6.embl.de/download/emapperdb-5.0.2/eggnog.taxa.tar.gz
    tar xzf eggnog.taxa.tar.gz -C eggnog_data
    """

    stub:
    """
    mkdir eggnog_data
    touch eggnog_data/eggnog.db
    touch eggnog_data/eggnog.taxa.db
    touch eggnog_data/eggnog.taxa.db.traverse.pkl
    touch eggnog_proteins.dmnd
    """
}
