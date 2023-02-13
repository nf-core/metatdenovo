process CAT {
    tag "${meta.id}-${db_name}"

    conda "bioconda::cat=4.6 bioconda::diamond=2.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cat:5.2.3--hdfd78af_1' :
        'my_image:latest' }"

    input:
    tuple val(meta), path(assembly)
    tuple val(db_name), path("database/*"), path("taxonomy/*")

    output:
    path("*.names.txt.gz")                    , emit: tax_classification
    path("raw/*.ORF2LCA.txt.gz")              , emit: orf2lca
    path("raw/*.predicted_proteins.faa.gz")   , emit: faa
    path("raw/*.predicted_proteins.gff.gz")   , emit: gff
    path("raw/*.log")                         , emit: log
    path("raw/*.contig2classification.txt.gz"), emit: tax_classification_taxids
    path "versions.yml"                       , emit: versions

    script:
    def official_taxonomy = params.cat_official_taxonomy ? "--only_official" : ""
    """
    CAT contigs -c "$assembly" -d database/ -t taxonomy/ -n 4 -o "${meta.id}"
    CAT add_names -i "${meta.id}.ORF2LCA.txt" -o "${meta.id}.ORF2LCA.names.txt" -t taxonomy/ ${official_taxonomy}
    CAT add_names -i "${meta.id}.contig2classification.txt" -o "${meta.id}.contig2classification.names.txt" -t taxonomy/ ${official_taxonomy}

    mkdir raw
    mv *.ORF2LCA.txt *.predicted_proteins.faa *.predicted_proteins.gff *.log *.bin2classification.txt raw/
    gzip \
        "raw/${meta.id}.concatenated.predicted_proteins.faa" \
        "raw/${meta.id}.concatenated.predicted_proteins.gff" \
        "raw/${meta.id}.bin2classification.txt" \
        "${meta.id}.ORF2LCA.names.txt" \
        "${meta.id}.bin2classification.names.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        CAT: \$(CAT --version | sed "s/CAT v//; s/(.*//")
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
