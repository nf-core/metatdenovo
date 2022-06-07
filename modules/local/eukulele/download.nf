process EUKULELE_DB {
    tag '$meta.id'
    label 'process_long'
    
    conda (params.enable_conda ? "bioconda::eukulele=2.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:2.0.3--pyh723bec7_0' :
        'quay.io/biocontainers/eukulele:2.0.3--pyh723bec7_0' }"

    input:
    
    output:
    path "versions.yml", emit: version
    path("phylodb/")   , emit: phylo_db , optional: true
    path("mmetsp/")    , emit: mmetsp_db, optional: true

    
    script:
    def args = task.ext.args ?: ''

    """
    EUKulele \\
    download \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        eukulele: \$(echo \$(EUKulele --version 2>&1) | sed 's/Running EUKulele with command line arguments, as no valid configuration file was provided.//; s/The current EUKulele version is//g')
    END_VERSIONS

    """
}
