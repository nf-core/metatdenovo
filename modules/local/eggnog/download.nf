process EGGNOG_DOWNLOAD {
    tag '$meta.id'
    label 'process_low'

    conda "bioconda::eggnog-mapper=2.1.9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.9--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.9--pyhdfd78af_0' }"


    output:
    path("eggnog")     , emit: db
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    """

    mkdir eggnog

    download_eggnog_data.py \\
        $args \\
        -y \\
        --data_dir eggnog

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
