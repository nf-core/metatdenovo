//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_SCAN     } from '../../modules/local/kofamscan/scan'
include { KOFAMSCAN_DOWNLOAD } from '../../modules/local/kofamscan/download'
include { SUM_KOFAMSCAN      } from '../../modules/local/sum_kofamscan'

workflow KOFAMSCAN {

    take:
    kofamscan // Channel: val(meta), path(fasta)
    fcs       // featureCounts output

    main:
    ch_versions = Channel.empty()

    KOFAMSCAN_DOWNLOAD()
    ch_versions = ch_versions.mix ( KOFAMSCAN_DOWNLOAD.out.versions )

    KOFAMSCAN_SCAN( kofamscan, KOFAMSCAN_DOWNLOAD.out.ko_list, KOFAMSCAN_DOWNLOAD.out.koprofiles )
    ch_versions = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)

    SUM_KOFAMSCAN( KOFAMSCAN_SCAN.out.kout, fcs )
    ch_versions = ch_versions.mix(SUM_KOFAMSCAN.out.versions)

    emit:
    kofam_table_out   = KOFAMSCAN_SCAN.out.kout
    kofamscan_summary = SUM_KOFAMSCAN.out.kofamscan_summary
    versions          = ch_versions
}
