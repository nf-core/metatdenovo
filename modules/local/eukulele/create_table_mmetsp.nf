process TABLE_MMETSP {
    tag '$meta.id'
    label 'process_low'
    conda (params.enable_conda ? "bioconda::eukulele=1.0.6-0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/eukulele:1.0.6--pyh723bec7_0' :
        'quay.io/biocontainers/eukulele:1.0.6--pyh723bec7_0' }"

    input:
    path(db)
    path(taxatable)

    output:
    path("./eukulele/"), emit: db  

    script:
    def args = task.ext.args ?: ''
    target  = params.eukulele_dbpath + "/" + params.eukulele_db + "/prot-map.json"
    target2 = params.eukulele_dbpath + "/" + params.eukulele_db + "/tax-table.txt"
    target3 = params.eukulele_dbpath + "/" + params.eukulele_db + "/reference.pep.fa"
    target4 = params.eukulele_dbpath + "/" + params.eukulele_db + "/taxonomy-table.txt"
    db = params.eukulele_db

    """
    if [ ! -e $target ]; then
        create_protein_table.py  \\
            --infile_peptide $db \\
            --infile_taxonomy $taxatable \\
            --outfile_json $target \\
            --output $target2 \\
            --delim "\t" \\
            --col_source_id strain_name \\
            --taxonomy_col_id taxonomy \\
            --column 2 \\
            --reformat_tax
        echo "$target and $target2 are now available" >> eukulele_create_table.log
    else
        echo "$target and $target2 already present" >> eukulele_create_table.log
    fi

    mkdir -p eukulele/$db
    ln -s $target ./eukulele/$db
    ln -s $target2 ./eukulele/$db
    ln -s $target3 ./eukulele/$db/
    ln -s $target4 ./eukulele/$db/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       eukulele: \$(echo \$(EUKulele --version 2>&1) | sed 's/Running EUKulele with command line arguments, as no valid configuration file was provided.//; s/The current EUKulele version is//g')
    END_VERSIONS

    """
}
