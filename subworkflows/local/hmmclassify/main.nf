include { HMMER_HMMSEARCH  } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { HMMER_HMMRANK    } from '../../../modules/nf-core/hmmer/hmmrank/main'
include { SEQTK_HMMHITFAAS } from '../../../modules/local/seqtk/hmmhitfaas/main'

workflow HMMCLASSIFY {

    take:
    ch_hmmclassify // channel: [ val(meta), [ hmm, aa_fasta ] ]

    main:

    HMMER_HMMSEARCH (
        ch_hmmclassify
            .map { meta, hmm, seqdb -> [ [ id: "${meta.id}.${hmm.baseName}", caller: meta.caller ], hmm, seqdb, false, true, false ] }
    )

    // Group hmm-search summaries by caller before ranking, so a caller's HMMER_HMMRANK invocation
    // only ever sees its own hmm files' summaries -- not another simultaneously-active caller's.
    HMMER_HMMRANK (
        HMMER_HMMSEARCH.out.target_summary
            .map { meta, summary -> [ meta.caller, summary ] }
            .groupTuple()
            .join(
                ch_hmmclassify
                    .map { meta, _hmm, _seqdb -> [ meta.caller, meta ] }
                    .distinct()
            )
            .map { _caller, summaries, callerMeta -> [ callerMeta, summaries ] }
    )

    SEQTK_HMMHITFAAS(
        HMMER_HMMRANK.out.hmmrank
            .join(ch_hmmclassify)
            .map { meta, hmmrank, _hmms, faa -> [ meta, hmmrank, faa ] }
    )

    emit:
    hmmrank  = HMMER_HMMRANK.out.hmmrank
    faas     = SEQTK_HMMHITFAAS.out.faas
}
