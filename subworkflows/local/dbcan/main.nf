//
// Run dbCAN CAZyme annotation on called ORFs, first optionally downloading the required database
//

include { RUNDBCAN_DATABASE         } from '../../../modules/nf-core/rundbcan/database/main'
include { RUNDBCAN_CAZYMEANNOTATION } from '../../../modules/nf-core/rundbcan/cazymeannotation/main'
include { DBCAN_FORMAT              } from '../../../modules/local/dbcan/format/main'
include { DBCAN_SUM                 } from '../../../modules/local/dbcan/sum/main'

workflow DBCAN {
    take:
    faa
    collect_fcs

    main:

    RUNDBCAN_DATABASE()

    RUNDBCAN_CAZYMEANNOTATION(faa, RUNDBCAN_DATABASE.out.dbcan_db)

    DBCAN_FORMAT(RUNDBCAN_CAZYMEANNOTATION.out.cazyme_annotation)

    DBCAN_SUM(DBCAN_FORMAT.out.dbcantsv, collect_fcs)

    emit:
    cazyme_annotation = DBCAN_FORMAT.out.dbcantsv
    sumtable          = DBCAN_SUM.out.dbcan_summary
}
