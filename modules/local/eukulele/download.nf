process WGET_DB {
    tag '$meta.id'
    label 'process_low'
    
    input:
    val(fasta_url)
    val(tax_url)

    output:
    path("./*")        , emit: database
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    
    target = params.eukulele_dbpath + "/" + params.eukulele_db + "/reference.pep.fa.gz"
    //file(params.eukulele_paths[params.eukulele_db]['fasta_url']) ? target : ''
    //file(params.eukulele_paths[params.eukulele_db]['fasta_url'].endWith('.gz')) ? target | ".gz" : target
    //fetchcmd = "wget -O $target $fasta_url" && fasta_url.endWith('.gz') ? target | "; unpigz $target" : target
    fetchcmd = "wget -O $target $fasta_url"
    
    target2 = params.eukulele_dbpath + "/" + params.eukulele_db + "/taxonomy-table.txt"
    //file(params.eukulele_paths[params.eukulele_db]['tax_url']) ? target2 : ''
    //file(params.eukulele_paths[params.eukulele_db]['tax_url'].endWith('.gz')) ? target2 | ".gz" : target2
    //getchcmd = "wget -O $target2 $tax_url" && tax_url.endWith('.gz') ? target | "; unpigz $target2" : target2
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
