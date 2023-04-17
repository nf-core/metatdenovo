//
// Run prodigal as orf caller then generate nice format for gff
//

include { PRODIGAL as PRODIGAL_MODULE } from '../../modules/nf-core/prodigal/main'
include { FORMAT_PRODIGAL_GFF         } from '../../modules/local/format_prodigal_gff.nf'
include { FORMAT_PRODIGAL_FAA         } from '../../modules/local/format_prodigal_faa.nf'

workflow PRODIGAL {
    take:
        fastafile

    main:
        ch_versions = Channel.empty()
        
        PRODIGAL_MODULE     ( fastafile, 'gff' )
        FORMAT_PRODIGAL_GFF ( PRODIGAL_MODULE.out.gene_annotations )
        FORMAT_PRODIGAL_FAA ( PRODIGAL_MODULE.out.amino_acid_fasta )

        ch_versions = ch_versions.mix(PRODIGAL_MODULE.out.versions)

    emit:
        faa     = FORMAT_PRODIGAL_FAA.out.format_faa
        gff     = FORMAT_PRODIGAL_GFF.out.format_gff
        versions = ch_versions

}
