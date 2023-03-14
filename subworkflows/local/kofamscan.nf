//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_SCAN     } from '../../modules/local/kofamscan/scan'
include { KOFAMSCAN_DOWNLOAD } from '../../modules/local/kofamscan/download'

workflow KOFAMSCAN {

    take:
        kofamscan // Channel: val(meta), path(fasta)
        ko_list_file
        koprofiles_dir

    main:
        ch_versions = Channel.empty()
        
        KOFAMSCAN_DOWNLOAD ( ko_list_file, koprofiles_dir )
        ch_versions    = ch_versions.mix(KOFAMSCAN_DOWNLOAD.out.versions)

        KOFAMSCAN_SCAN( kofamscan, KOFAMSCAN_DOWNLOAD.out.ko_list, KOFAMSCAN_DOWNLOAD.out.ko_profiles)
        ch_versions    = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)

    emit:
        kofam_table_out = KOFAMSCAN_SCAN.out.kout
        versions        = ch_versions 
}
