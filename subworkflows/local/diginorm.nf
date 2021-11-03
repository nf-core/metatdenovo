//
// Perform digital normalization, i.e. subsetting of reads, with the khmer package.
//

params.diginorm_normalizebymedian       = [:]

include { KHMER_NORMALIZEBYMEDIAN } from '../../modules/local/khmer/normalizebymedian' addParams( options: params.diginorm_normalizebymedian )

workflow DIGINORM {
    take:
    pe_reads
    se_reads
    name

    main:
    ch_versions = Channel.empty()

    KHMER_NORMALIZEBYMEDIAN(pe_reads, se_reads, name)
    ch_versions = ch_versions.mix(KHMER_NORMALIZEBYMEDIAN.out.versions)

    emit:
    reads    = KHMER_NORMALIZEBYMEDIAN.out.reads
    versions = ch_versions
}
