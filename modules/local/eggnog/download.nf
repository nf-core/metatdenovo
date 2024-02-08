process EGGNOG_DOWNLOAD {
    tag 'EggNOG'
    label 'process_low'

    conda "bioconda::eggnog-mapper=2.1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.9--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0' }"

    input:
    //path "eggnog_dbpath"

    output:
    //path("eggnog_db")  , emit: db
    path('./.'), emit: db
    path('*'), emit: files
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """

    #mkdir eggnog_db

    download_eggnog_data.py \\
        $args \\
        -y \\
        --data_dir .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//' | sed 's/ \\/ Expected.*//')
    END_VERSIONS

    """

    stub:

    """

    mkdir eggnog
    touch ./eggnog/eggnog.db
    touch ./eggnog/eggnog.taxa.db
    touch ./eggnog/eggnog.taxa.db.traverse.pkl
    ln -s eggnog/* ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//' )
    END_VERSIONS

    """
}
