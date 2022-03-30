process WGET_DB {
    tag '$meta.id'
    label 'process_long'
    
    conda (params.enable_conda ? "bioconda::eukulele=1.0.6-0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:1.0.6--pyh723bec7_0' :
        'quay.io/biocontainers/eukulele:1.0.6--pyh723bec7_0' }"

    input:
    val(fastafile)
    val(taxonomyfile)
    
    output:
    path "versions.yml"         , emit: version
    path("./reference.pep.fa")  , emit: database
    path("./taxonomy-table.txt"), emit: taxonomy
    
    script:
    def args = task.ext.args ?: ''
    target = params.eukulele_dbpath + "/" + params.eukulele_db + "/reference.pep.fa"
    target1 = params.eukulele_dbpath + "/" + params.eukulele_db + "/reference.pep.fa.gz"
    target2 = params.eukulele_dbpath + "/" + params.eukulele_db + "/taxonomy-table.txt"
    target3 = params.eukulele_dbpath + "/" + params.eukulele_db + "/taxonomy-table.txt.gz"
    db = params.eukulele_db 

    """
    if [ ! -e $target ]; then
        wget -O $target1 $fastafile
        gunzip $target1
        echo "Downloaded $fastafile to $target" >> eukulele_dbdwnl.log
    else
        echo "$target already present" >> eukulele_dbdwnl.log
    fi
    
    if [ ! -e $target2 ]; then
        wget -O $target3 $taxonomyfile
        gunzip $target3
        echo "Downloaded $taxonomyfile  to $target2" >> eukulele_dbdwnl.log
    else
        echo "$target2 already present" >> eukulele_dbdwnl.log
    fi
    
    ln -s $target ./
    ln -s $target2 ./
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo \$(wget --version 2>&1))
    END_VERSIONS
    """
}
