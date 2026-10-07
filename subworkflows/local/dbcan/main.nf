//
// Run dbCAN CAZyme annotation on called ORFs, first optionally downloading the required database
//

include { RUNDBCAN_DATABASE         } from '../../../modules/nf-core/rundbcan/database/main'
include { RUNDBCAN_CAZYMEANNOTATION } from '../../../modules/nf-core/rundbcan/cazymeannotation/main'
include { DBCAN_FORMAT              } from '../../../modules/local/dbcan/format/main'
include { DBCAN_SUM                 } from '../../../modules/local/dbcan/sum/main'

workflow DBCAN {
    take:
    faa            // channel: [ val(meta), path(faa) ]
    feature_counts // channel: [ val(meta), path(fcs) ] -- meta.caller must match faa's

    main:

    RUNDBCAN_DATABASE()

    RUNDBCAN_CAZYMEANNOTATION(faa, RUNDBCAN_DATABASE.out.dbcan_db)

    DBCAN_FORMAT(RUNDBCAN_CAZYMEANNOTATION.out.cazyme_annotation)

    ch_dbcan_sum_input = DBCAN_FORMAT.out.dbcantsv
        .map { meta, tsv -> [ meta.caller, meta, tsv ] }
        .join( feature_counts.map { meta, fcs -> [ meta.caller, fcs ] } )
        .map { _caller, meta, tsv, fcs -> [ meta, tsv, fcs ] }

    DBCAN_SUM(ch_dbcan_sum_input)

    emit:
    cazyme_annotation = DBCAN_FORMAT.out.dbcantsv
    sumtable          = DBCAN_SUM.out.dbcan_summary
}
