//
// TRANSDECODER SUBWORKFLOW
//

include { TRANSDECODER_LONGORF as LONGORF } from '../../modules/nf-core/transdecoder/longorf/main'
include { TRANSDECODER_PREDICT as PREDICT } from '../../modules/nf-core/transdecoder/predict/main'

workflow TRANSDECODER {
    take:
    contigs       // channel: [ val(meta), [ contigs ] ]

    main:
    ch_versions  = Channel.empty()

    LONGORF (contigs)
    ch_versions = ch_versions.mix(LONGORF.out.versions)
    PREDICT (contigs, LONGORF.out.folder)
    ch_versions = ch_versions.mix(PREDICT.out.versions)

    emit:
    gff      = PREDICT.out.gff3
    cds      = PREDICT.out.cds
    pep      = PREDICT.out.pep
    bed      = PREDICT.out.bed

    versions = ch_versions

}
