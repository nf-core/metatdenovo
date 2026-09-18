//
// Run InterProScan on called ORFs, first downloading the required database
//

include { INTERPROSCAN_DATABASE      } from '../../../modules/local/interproscan/download/main'
include { INTERPROSCAN_STRIPASTERISK } from '../../../modules/local/interproscan/stripasterisk/main'
include { INTERPROSCAN               } from '../../../modules/nf-core/interproscan/main'
include { INTERPROSCAN_FORMAT        } from '../../../modules/local/interproscan/format/main'
include { INTERPROSCAN_SUM           } from '../../../modules/local/interproscan/sum/main'

workflow INTERPROSCAN_ANNOTATION {
    take:
    faa            // channel: [ val(meta), path(faa) ]
    feature_counts // channel: [ val(meta), path(fcs) ] -- meta.caller must match faa's

    main:

    INTERPROSCAN_DATABASE(params.interproscan_db_url)

    INTERPROSCAN_STRIPASTERISK(faa)

    INTERPROSCAN(INTERPROSCAN_STRIPASTERISK.out.faa, INTERPROSCAN_DATABASE.out.db)

    INTERPROSCAN_FORMAT(INTERPROSCAN.out.tsv)

    ch_interproscan_sum_input = INTERPROSCAN_FORMAT.out.tsv
        .map { meta, tsv -> [ meta.caller, meta, tsv ] }
        .join( feature_counts.map { meta, fcs -> [ meta.caller, fcs ] } )
        .map { _caller, meta, tsv, fcs -> [ meta, tsv, fcs ] }

    INTERPROSCAN_SUM(ch_interproscan_sum_input)

    emit:
    annotation = INTERPROSCAN_FORMAT.out.tsv
    sumtable   = INTERPROSCAN_SUM.out.interproscan_summary
}
