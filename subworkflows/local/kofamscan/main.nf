//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_DOWNLOAD } from '../../../modules/local/kofamscan/download/main'
include { KOFAMSCAN_SCAN     } from '../../../modules/local/kofamscan/scan/main'
include { KOFAMSCAN_UNIQUE   } from '../../../modules/local/kofamscan/unique/main'
include { KOFAMSCAN_SUM      } from '../../../modules/local/kofamscan/sum/main'

workflow KOFAMSCAN {

    take:
    kofamscan       // Channel: val(meta), path(fasta)
    fcs             // featureCounts output
    ko_list_url     // string: URL to download the KOfam ko_list file from
    profiles_url    // string: URL to download the KOfam HMM profiles archive from

    main:

    KOFAMSCAN_DOWNLOAD( ko_list_url, profiles_url )

    KOFAMSCAN_SCAN(kofamscan, KOFAMSCAN_DOWNLOAD.out.ko_list, KOFAMSCAN_DOWNLOAD.out.koprofiles)

    KOFAMSCAN_UNIQUE(KOFAMSCAN_SCAN.out.kofamtsv)

    KOFAMSCAN_SUM(KOFAMSCAN_SCAN.out.kout, fcs)

    emit:
    kofam_table_out   = KOFAMSCAN_SCAN.out.kout
    kofam_table_tsv   = KOFAMSCAN_SCAN.out.kofamtsv
    kofam_table_uniq  = KOFAMSCAN_UNIQUE.out.kofamuniq
    kofamscan_summary = KOFAMSCAN_SUM.out.kofamscan_summary
}
