//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_SCAN      } from '../../modules/local/kofamscan/scan'
include { DOWNLOAD_KOFMASCAN_DB as DOWNLOAD } from '../../modules/local/kofamscan/download'

workflow KOFAMSCAN {

    take:
        kofamscan // Channel: val(meta), path(fasta)
        databases // Channel: path(ko_list), path(koprofiles)
        check_db  // Channel: path(ko_list)

    main:
        ch_versions = Channel.empty()
        
        if ( ! check_db.exists() ) {
            DOWNLOAD ( )
            ch_ko_profiles = DOWNLOAD.out.profiles
            ch_ko_list     = DOWNLOAD.out.ko_list
            ch_versions    = ch_versions.mix(DOWNLOAD.out.versions)
            ch_ko_db = ch_ko_list
                .map { [ [id: 'ko_database'], it ] }
                .combine( ch_ko_profiles)
            KOFAMSCAN_SCAN( kofamscan, ch_ko_db )
            ch_versions = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)
        } else {
            KOFAMSCAN_SCAN( kofamscan, ch_ko_list, ch_ko_profiles )
            ch_versions = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)
        }

    emit:
        kofam_table_out = KOFAMSCAN_SCAN.out.kout
        versions        = ch_versions 
}
