//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_SCAN     } from '../../modules/local/kofamscan/scan'
include { KOFAMSCAN_DOWNLOAD } from '../../modules/local/kofamscan/download'
include { SUM_KOFAMSCAN      } from '../../modules/local/sum_kofamscan'

workflow KOFAMSCAN {

    take:
        kofamscan // Channel: val(meta), path(fasta)
        kofam_dir
        fcs

    main:
        ch_versions = Channel.empty()
        
        KOFAMSCAN_DOWNLOAD ( kofam_dir )

        KOFAMSCAN_SCAN( kofamscan, KOFAMSCAN_DOWNLOAD.out.ko_list, KOFAMSCAN_DOWNLOAD.out.koprofiles )
        ch_versions = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)
        
        SUM_KOFAMSCAN( KOFAMSCAN_SCAN.out.kout, fcs )


    emit:
        kofam_table_out   = KOFAMSCAN_SCAN.out.kout
        kofamscan_summary = SUM_KOFAMSCAN.out.kofamscan_summary
        versions          = ch_versions
}
