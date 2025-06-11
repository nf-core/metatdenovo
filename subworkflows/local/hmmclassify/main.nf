include { HMMER_HMMSEARCH  } from '../../../modules/nf-core/hmmer/hmmsearch/main'
include { HMMRANK          } from '../../../modules/local/hmmrank/main'
include { SEQTK_HMMHITFAAS } from '../../../modules/local/seqtk/hmmhitfaas/main'

workflow HMMCLASSIFY {

    take:
    ch_hmmclassify // channel: [ val(meta), [ hmm, aa_fasta ] ]

    main:
    ch_versions = Channel.empty()

    HMMER_HMMSEARCH (
        ch_hmmclassify
            .map { meta, hmm, seqdb -> [ [ id: "${meta.id}.${hmm.baseName}" ], hmm, seqdb, false, true, false ] }
    )
    ch_versions = ch_versions.mix(HMMER_HMMSEARCH.out.versions.first())

    HMMRANK (
        ch_hmmclassify
            .map { meta, hmm, seqdb -> meta }
            .distinct()
            .combine ( HMMER_HMMSEARCH.out.target_summary.collect { meta, summary -> summary } )
            .map { [ it[0], it[1..-1] ] }
    )
    ch_versions = ch_versions.mix(HMMRANK.out.versions)

    SEQTK_HMMHITFAAS(
        HMMRANK.out.hmmrank
            .join(ch_hmmclassify)
            .map { meta, hmmrank, hmms, faa -> [ meta, hmmrank, faa ] }
    )
    ch_versions = ch_versions.mix(SEQTK_HMMHITFAAS.out.versions)

    emit:
    HMMRANK.out.hmmrank
    SEQTK_HMMHITFAAS.out.faas

    versions = ch_versions                     // channel: [ versions.yml ]
}
