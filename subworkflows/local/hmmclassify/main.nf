include { HMMER_HMMSEARCH  } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { HMMRANK          } from '../../../modules/local/hmmrank/main'
include { SEQTK_HMMHITFAAS } from '../../../modules/local/seqtk/hmmhitfaas/main'

workflow HMMCLASSIFY {

    take:
    ch_hmmclassify // channel: [ val(meta), [ hmm, aa_fasta ] ]

    main:

    HMMER_HMMSEARCH (
        ch_hmmclassify
            .map { meta, hmm, seqdb -> [ [ id: "${meta.id}.${hmm.baseName}" ], hmm, seqdb, false, true, false ] }
    )

    HMMRANK (
        ch_hmmclassify
            .map { meta, _hmm, _seqdb -> meta }
            .distinct()
            .combine ( HMMER_HMMSEARCH.out.target_summary.collect { _meta, summary -> summary } )
            .map { summary -> [ summary[0], summary[1..-1] ] }
    )

    SEQTK_HMMHITFAAS(
        HMMRANK.out.hmmrank
            .join(ch_hmmclassify)
            .map { meta, hmmrank, _hmms, faa -> [ meta, hmmrank, faa ] }
    )

    emit:
    hmmrank  = HMMRANK.out.hmmrank
    faas     = SEQTK_HMMHITFAAS.out.faas
}
