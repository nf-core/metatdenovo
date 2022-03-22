process WGET_DB {
    tag '$meta.id'
    label 'process_low'
    
    input:
    val(fasta_url)
    val(tax_url)

    output:
    path("~/eukulele") , emit: database
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    
    target = params.eukulele_dbpath + "/" + params.eukulele_db + "/reference.pep.fa.gz"
    fetchcmd = "wget -O $target $fasta_url"
    
    target2 = params.eukulele_dbpath + "/" + params.eukulele_db + "/taxonomy-table.txt"
    getchcmd = "wget -O $target2 $tax_url"
    
    """ 
    $fetchcmd
    $getchcmd


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(echo \$(wget --version 2>&1))
    END_VERSIONS
    """
}
