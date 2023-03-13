//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN as EXEC_ANNOTATION      } from '../../modules/local/kofamscan/main'
include { DOWNLOAD_KOFMASCAN_DB as DOWNLOAD } from '../../modules/local/kofamscan/download'

workflow KOFAMSCAN {

    take:
        kofamscan // Channel: val(meta), path(fasta)
        databases // Channel: val(meta), path(ko_list), path(koprofiles)
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
            EXEC_ANNOTATION( kofamscan, ch_ko_db )
            ch_versions = ch_versions.mix(EXEC_ANNOTATION.out.versions)
        } else {
            EXEC_ANNOTATION( kofamscan, databases )
            ch_versions = ch_versions.mix(EXEC_ANNOTATION.out.versions)
        }

    emit:
        kofam_table_out = EXEC_ANNOTATION.out.kout
        versions        = ch_versions 
}
