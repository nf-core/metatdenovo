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

    // TransDecoder.LongOrfs works per-transcript, so splitting into batches (same
    // rationale/mechanism as PROKKA_SUBSETS/METAEUK) is safe there. TransDecoder.Predict is a
    // different story: it self-trains a coding/noncoding hexamer model from the top longest
    // ORFs in whatever input it's given (-T, default 500 of the top 5000 candidates), and its
    // CLI has no option to reuse a model trained elsewhere. Batching therefore means each
    // batch trains its own model from its own slice of the assembly rather than one model
    // over the whole thing -- the same ORF sequence can in principle score slightly
    // differently depending on which batch it lands in. A large --transdecoder_batchsize
    // (bigger than the assembly) recovers the unbatched, single-model behaviour exactly; the
    // default trades some of that consistency for the same resilience/resume benefits
    // Prokka and MetaEuk batching already get -- a killed/failed batch only costs re-running
    // that batch, not a multi-hour whole-assembly job.
    ch_batches = fasta
        .map { _meta, contigs -> contigs }
        .splitFasta(size: batchsize, file: true)
        .map { ctg -> [ [ id: ctg.getBaseName() ], ctg ] }

    TRANSDECODER_LONGORF ( ch_batches )

    // TRANSDECODER_LONGORF.out.folder carries no meta, so pairing it with ch_batches
    // positionally (as two separate process inputs) relies on both channels emitting in the
    // same order -- which isn't guaranteed once batches run concurrently: a smaller/faster
    // batch's LongOrfs can finish before an earlier, bigger batch's, silently swapping which
    // fasta gets paired with which batch's LongOrfs directory. Join explicitly on id instead,
    // recovered from the folder's own parent directory name (LONGORF names it
    // "${meta.id}/${fasta_no_gz}.transdecoder_dir").
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
