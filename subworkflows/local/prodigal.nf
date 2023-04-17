//
// Run prodigal as orf caller then generate nice format for gff
//

include { PRODIGAL as PRODIGAL_MODULE } from '../../modules/nf-core/prodigal/main'
include { FORMAT_PRODIGAL_GFF         } from '../../modules/local/format_prodigal'

workflow PRODIGAL {
    take:
        fastafile

    main:
        ch_versions = Channel.empty()

        PRODIGAL_MODULE     ( fastafile, 'gff' )
        FORMAT_PRODIGAL_GFF ( PRODIGAL_MODULE.out.gene_annotations )
        ch_versions = ch_versions.mix(PRODIGAL_MODULE.out.versions)

    emit:
        faa     = PRODIGAL_MODULE.out.amino_acid_fasta
        gff     = FORMAT_PRODIGAL_GFF.out.format_gff
        versions = ch_versions

}
