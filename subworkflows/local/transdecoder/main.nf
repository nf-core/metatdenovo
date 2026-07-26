//
// Call ORFs with TransDecoder: TransDecoder.LongOrfs then TransDecoder.Predict
//

include { TRANSDECODER_LONGORF } from '../../../modules/nf-core/transdecoder/longorf/main'
include { TRANSDECODER_PREDICT } from '../../../modules/nf-core/transdecoder/predict/main'

workflow TRANSDECODER {

    take:
    fasta // channel: [ val(meta), path(fasta) ]

    main:

    TRANSDECODER_LONGORF ( fasta )

    TRANSDECODER_PREDICT ( fasta, TRANSDECODER_LONGORF.out.folder )

    emit:
    pep = TRANSDECODER_PREDICT.out.pep // channel: [ val(meta), path(pep)  ]
    gff = TRANSDECODER_PREDICT.out.gff3 // channel: [ val(meta), path(gff3) ]
    cds = TRANSDECODER_PREDICT.out.cds // channel: [ val(meta), path(cds)  ]
    bed = TRANSDECODER_PREDICT.out.bed // channel: [ val(meta), path(bed)  ]
}
