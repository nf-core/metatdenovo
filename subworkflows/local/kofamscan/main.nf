//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_DOWNLOAD }                        from '../../../modules/local/kofamscan/download/main'
include { PIGZ_UNCOMPRESS as UNPIGZ_KOFAMSCAN_FASTA } from '../../../modules/nf-core/pigz/uncompress/main'
include { KOFAMSCAN as KOFAMSCAN_SCAN }                from '../../../modules/nf-core/kofamscan/main'
include { PIGZ_COMPRESS as PIGZ_KOFAMSCAN_RAW }        from '../../../modules/nf-core/pigz/compress/main'
include { KOFAMSCAN_FORMAT }                          from '../../../modules/local/kofamscan/format/main'
include { KOFAMSCAN_UNIQUE }                          from '../../../modules/local/kofamscan/unique/main'
include { KOFAMSCAN_SUM }                             from '../../../modules/local/kofamscan/sum/main'

workflow KOFAMSCAN {

    take:
    kofamscan       // Channel: val(meta), path(fasta)
    fcs             // featureCounts output
    ko_list_url     // string: URL to download the KOfam ko_list file from
    profiles_url    // string: URL to download the KOfam HMM profiles archive from

    main:

    KOFAMSCAN_DOWNLOAD( ko_list_url, profiles_url )

    UNPIGZ_KOFAMSCAN_FASTA( kofamscan )

    KOFAMSCAN_SCAN( UNPIGZ_KOFAMSCAN_FASTA.out.file, KOFAMSCAN_DOWNLOAD.out.koprofiles, KOFAMSCAN_DOWNLOAD.out.ko_list )

    PIGZ_KOFAMSCAN_RAW( KOFAMSCAN_SCAN.out.tsv )

    KOFAMSCAN_FORMAT( KOFAMSCAN_SCAN.out.tsv )

    KOFAMSCAN_UNIQUE( KOFAMSCAN_FORMAT.out.kofamtsv )

    KOFAMSCAN_SUM( PIGZ_KOFAMSCAN_RAW.out.archive, fcs )

    emit:
    kofam_table_out   = PIGZ_KOFAMSCAN_RAW.out.archive
    kofam_table_tsv   = KOFAMSCAN_FORMAT.out.kofamtsv
    kofam_table_uniq  = KOFAMSCAN_UNIQUE.out.kofamuniq
    kofamscan_summary = KOFAMSCAN_SUM.out.kofamscan_summary
}
