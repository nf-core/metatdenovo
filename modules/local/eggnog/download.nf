process EGGNOG_DOWNLOAD {
    tag '$meta.id'
    label 'process_low'

    conda "bioconda::eggnog-mapper=2.1.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.6--pyhdfd78af_0':
        'biocontainers/eggnog-mapper:2.1.6--pyhdfd78af_0' }"

    input:
    path(dbpath)

    output:
    path("$dbpath/*.db")          , emit: eggnog_db
    path("$dbpath/*.taxa.db")     , emit: eggnog_taxa
    path("$dbpath/*.pkl")         , emit: eggnog_traverse
    path("$dbpath/")              , emit: db
    path("$dbpath/*.dmnd")        , emit: proteins, optional: true
    path("$dbpath/hmmer/")        , emit: hmmer   , optional: true
    path("$dbpath/mmseqs/")       , emit: mmseqs  , optional: true
    path("$dbpath/pfam/")         , emit: pfam    , optional: true
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """

    download_eggnog_data.py \\
        $args \\
        -y \\
        --data_dir $dbpath

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eggnog: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//' )
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
