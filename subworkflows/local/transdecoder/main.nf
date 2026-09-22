//
// Call ORFs with TransDecoder (LongOrfs then Predict) on transcripts split by size, then
// concatenate output
//

include { TRANSDECODER_LONGORF                     } from '../../../modules/nf-core/transdecoder/longorf/main'
include { TRANSDECODER_PREDICT                     } from '../../../modules/nf-core/transdecoder/predict/main'
include { FIND_CONCATENATE as TRANSDECODER_PEP_CAT } from '../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as TRANSDECODER_GFF_CAT } from '../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as TRANSDECODER_CDS_CAT } from '../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as TRANSDECODER_BED_CAT } from '../../../modules/nf-core/find/concatenate/main'

workflow TRANSDECODER {

    take:
    fasta     // channel: [ val(meta), path(fasta) ]
    batchsize // channel: strings like '10.MB'. Usually from params.transdecoder_batchsize

    main:

    // Predict self-trains per batch, so a batch size above the assembly size recovers unbatched results.
    ch_batches = fasta
        .map { _meta, contigs -> contigs }
        .splitFasta(size: batchsize, file: true)
        .map { ctg -> [ [ id: ctg.getBaseName() ], ctg ] }

    TRANSDECODER_LONGORF ( ch_batches )

    // LongOrfs output has no meta and emits out of order; join on the id from its parent dir name.
    ch_predict_in = ch_batches
        .map { meta, ctg -> [ meta.id, meta, ctg ] }
        .join( TRANSDECODER_LONGORF.out.folder.map { folder -> [ folder.getParent().getName(), folder ] } )
        .map { _id, meta, ctg, folder -> [ meta, ctg, folder ] }

    TRANSDECODER_PREDICT (
        ch_predict_in.map { meta, ctg, _folder -> [ meta, ctg ] },
        ch_predict_in.map { _meta, _ctg, folder -> folder }
    )

    ch_pep = fasta.map { meta, _contigs -> meta }
        .combine(TRANSDECODER_PREDICT.out.pep.collect { _meta, pep -> pep }.map { peps -> [ peps ] })
    TRANSDECODER_PEP_CAT(ch_pep)

    ch_gff = fasta.map { meta, _contigs -> meta }
        .combine(TRANSDECODER_PREDICT.out.gff3.collect { _meta, gff -> gff }.map { gffs -> [ gffs ] })
    TRANSDECODER_GFF_CAT(ch_gff)

    ch_cds = fasta.map { meta, _contigs -> meta }
        .combine(TRANSDECODER_PREDICT.out.cds.collect { _meta, cds -> cds }.map { cdss -> [ cdss ] })
    TRANSDECODER_CDS_CAT(ch_cds)

    ch_bed = fasta.map { meta, _contigs -> meta }
        .combine(TRANSDECODER_PREDICT.out.bed.collect { _meta, bed -> bed }.map { beds -> [ beds ] })
    TRANSDECODER_BED_CAT(ch_bed)

    emit:
    pep = TRANSDECODER_PEP_CAT.out.file_out // channel: [ val(meta), path(pep)  ]
    gff = TRANSDECODER_GFF_CAT.out.file_out // channel: [ val(meta), path(gff3) ]
    cds = TRANSDECODER_CDS_CAT.out.file_out // channel: [ val(meta), path(cds)  ]
    bed = TRANSDECODER_BED_CAT.out.file_out // channel: [ val(meta), path(bed)  ]
}
