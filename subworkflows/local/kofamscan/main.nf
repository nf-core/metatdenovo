//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_DOWNLOAD }         from '../../../modules/local/kofamscan/download/main'
include { KOFAMSCAN as KOFAMSCAN_SCAN } from '../../../modules/nf-core/kofamscan/main'
include { KOFAMSCAN_FORMAT }           from '../../../modules/local/kofamscan/format/main'
include { KOFAMSCAN_UNIQUE }           from '../../../modules/local/kofamscan/unique/main'
include { KOFAMSCAN_SUM }              from '../../../modules/local/kofamscan/sum/main'

workflow KOFAMSCAN {

    take:
    kofamscan       // Channel: val(meta), path(fasta)
    fcs             // channel: [ val(meta), path(fcs) ] -- meta.caller must match kofamscan's
    ko_list_url     // string: URL to download the KOfam ko_list file from
    profiles_url    // string: URL to download the KOfam HMM profiles archive from

    main:

    KOFAMSCAN_DOWNLOAD( ko_list_url, profiles_url )

    KOFAMSCAN_SCAN( kofamscan, KOFAMSCAN_DOWNLOAD.out.koprofiles, KOFAMSCAN_DOWNLOAD.out.ko_list )

    KOFAMSCAN_FORMAT( KOFAMSCAN_SCAN.out.tsv )

    KOFAMSCAN_UNIQUE( KOFAMSCAN_FORMAT.out.kofamtsv )

    ch_kofamscan_sum_input = KOFAMSCAN_SCAN.out.tsv
        .map { meta, tsv -> [ meta.caller, meta, tsv ] }
        .join( fcs.map { meta, f -> [ meta.caller, f ] } )
        .map { _caller, meta, tsv, f -> [ meta, tsv, f ] }

    KOFAMSCAN_SUM( ch_kofamscan_sum_input )

    emit:
    kofam_table_out   = KOFAMSCAN_SCAN.out.tsv
    kofam_table_tsv   = KOFAMSCAN_FORMAT.out.kofamtsv
    kofam_table_uniq  = KOFAMSCAN_UNIQUE.out.kofamuniq
    kofamscan_summary = KOFAMSCAN_SUM.out.kofamscan_summary
}
