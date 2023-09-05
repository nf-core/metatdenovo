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

        String directoryName = kofam_dir
        File directory       = new File(directoryName)
        String kofamdb       = directoryName + "ko_list"
        File kolistfile      = new File(kofamdb)

        if ( ! directory.exists() ) {
            directory.mkdir()
        }

        if ( ! kolistfile.exists() ) {
            KOFAMSCAN_DOWNLOAD ( kofam_dir )
            ch_dbpath   = KOFAMSCAN_DOWNLOAD.out.ko_list
            ch_profiles = KOFAMSCAN_DOWNLOAD.out.koprofiles

            ch_versions = ch_versions.mix ( KOFAMSCAN_DOWNLOAD.out.versions )
        } else {
            ch_dbpath   = Channel.fromPath(kolistfile)
            ch_profiles = Channel.fromPath(kofam_dir + "profiles")
        }

        KOFAMSCAN_SCAN( kofamscan, ch_dbpath, ch_profiles )
        ch_versions = ch_versions.mix(KOFAMSCAN_SCAN.out.versions)

        SUM_KOFAMSCAN( KOFAMSCAN_SCAN.out.kout, fcs )


    emit:
        kofam_table_out   = KOFAMSCAN_SCAN.out.kout
        kofamscan_summary = SUM_KOFAMSCAN.out.kofamscan_summary
        versions          = ch_versions
}
