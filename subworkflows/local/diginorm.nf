//
// Perform digital normalization, i.e. subsetting of reads, with the khmer package.
//

include { KHMER_NORMALIZEBYMEDIAN  } from '../../modules/local/khmer/normalizebymedian'
include { KHMER_FILTERABUND        } from '../../modules/local/khmer/filterabund'
include { KHMER_EXTRACTPAIREDREADS } from '../../modules/local/khmer/extractpairedreads'

workflow DIGINORM {
    take:
    pe_reads
    se_reads
    name

    main:
    ch_versions = Channel.empty()

    KHMER_NORMALIZEBYMEDIAN(pe_reads, se_reads, name)
    ch_versions = ch_versions.mix(KHMER_NORMALIZEBYMEDIAN.out.versions)

    KHMER_FILTERABUND(KHMER_NORMALIZEBYMEDIAN.out.reads, KHMER_NORMALIZEBYMEDIAN.out.graph, "${name}.nm")
    ch_versions = ch_versions.mix(KHMER_FILTERABUND.out.versions)

    KHMER_EXTRACTPAIREDREADS(KHMER_FILTERABUND.out.reads, "${name}.nm.fa")
    ch_versions = ch_versions.mix(KHMER_EXTRACTPAIREDREADS.out.versions)

    emit:
    pairs    = KHMER_EXTRACTPAIREDREADS.out.pairs
    singles  = KHMER_EXTRACTPAIREDREADS.out.singles
    versions = ch_versions
}
