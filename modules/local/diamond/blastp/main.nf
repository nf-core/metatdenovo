process DIAMOND_BLASTP {
    tag "${meta.id}.${meta2.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/diamond:2.1.10--h5ca1c30_3':
        'biocontainers/diamond:2.1.10--h5ca1c30_3' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(db)
    val outfmt
    val blast_columns

    output:
    tuple val(outmeta), path('*.blast*'), optional: true, emit: blast
    tuple val(outmeta), path('*.xml*')  , optional: true, emit: xml
    tuple val(outmeta), path('*.txt*')  , optional: true, emit: txt
    tuple val(outmeta), path('*.daa')   , optional: true, emit: daa
    tuple val(outmeta), path('*.sam*')  , optional: true, emit: sam
    tuple val(outmeta), path('*.tsv*')  , optional: true, emit: tsv
    tuple val(outmeta), path('*.paf*')  , optional: true, emit: paf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    outmeta = meta + [ db: meta2.id ]

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.${meta2.id}"
    def columns = blast_columns ? "${blast_columns}" : ''
    switch ( outfmt ) {
        case 0:   out_ext = "blast"; break
        case 5:   out_ext = "xml";   break
        case 6:   out_ext = "txt";   break
        case 100: out_ext = "daa";   break
        case 101: out_ext = "sam";   break
        case 102: out_ext = "tsv";   break
        case 103: out_ext = "paf";   break
        default:
            outfmt = 6;
            out_ext = 'txt';
            log.warn("Unknown output file format provided (${out_ext}): selecting DIAMOND default of tabular BLAST output (txt)");
            break
    }
    if ( args =~ /--compress\s+1/ ) out_ext += '.gz'
    """
    diamond \\
        blastp \\
        --threads ${task.cpus} \\
        --db ${db} \\
        --query ${fasta} \\
        --outfmt ${outfmt} ${columns} \\
        ${args} \\
        --out ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
