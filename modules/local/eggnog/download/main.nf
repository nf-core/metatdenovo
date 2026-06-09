process EGGNOG_DOWNLOAD {
    tag 'EggNOG'
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.9--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0' }"

    output:
    path "eggnog.db"                  , emit: eggnog_db
    path "eggnog_proteins.dmnd"       , emit: dmnd
    path "eggnog.taxa.db"             , emit: taxa_db
    path "eggnog.taxa.db.traverse.pkl", emit: pkl
    path "*"                          , emit: all
    tuple val("${task.process}"), val("eggnog-mapper"), eval('emapper.py --version | sed "s/.* emapper-//" | sed "s/ \\/ Expected.*//"'), emit: versions_emapper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    download_eggnog_data.py \\
        $args \\
        -y \\
        --data_dir .
    """

    stub:
    """
    mkdir eggnog
    touch ./eggnog/eggnog.db
    touch ./eggnog/eggnog.taxa.db
    touch ./eggnog/eggnog.taxa.db.traverse.pkl
    ln -s eggnog/* ./
    """
}
