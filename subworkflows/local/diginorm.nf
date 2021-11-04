//
// Perform digital normalization, i.e. subsetting of reads, with the khmer package.
//

params.diginorm_normalizebymedian_options       = [:]
params.diginorm_filterabund_options             = [:]

include { KHMER_NORMALIZEBYMEDIAN } from '../../modules/local/khmer/normalizebymedian' addParams( options: params.diginorm_normalizebymedian_options )
include { KHMER_FILTERABUND       } from '../../modules/local/khmer/filterabund'       addParams( options: params.diginorm_filterabund_options )

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

    emit:
    reads    = KHMER_FILTERABUND.out.reads
    versions = ch_versions
}
