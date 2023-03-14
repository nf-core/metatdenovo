//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_SCAN      } from '../../modules/local/kofamscan/scan'
include { KOFAMSCAN_DOWNLOAD } from '../../modules/local/kofamscan/download'

workflow KOFAMSCAN {

    take:
        kofamscan // Channel: val(meta), path(fasta)
        ko_list_file
        koprofiles_dir

    main:
        ch_versions = Channel.empty()
        
        // We're assuming that if the ko_list_file exists, the koprofiles directory is also downloaded.
        //if ( ! ko_list_file.exists() || ! koprofiles.HMMFILE.exists ) {
        //    KOFAMSCAN_DOWNLOAD ( )
        //    ch_versions    = ch_versions.mix(KOFAMSCAN_DOWNLOAD.out.versions)
        //    ch_ko_db = ch_ko_list
        //        .map { [ [id: 'ko_database'], it ] }
        //        .combine( ch_ko_profiles)
        //    KOFAMSCAN_SCAN( kofamscan, ch_ko_db )
        //    ch_versions = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)
        //} else {
        KOFAMSCAN_SCAN( kofamscan, ko_list_file, koprofiles_dir )
        ch_versions = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)
        //}

    emit:
        kofam_table_out = KOFAMSCAN_SCAN.out.kout
        versions        = ch_versions 
}
