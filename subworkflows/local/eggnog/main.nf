//
// Run eggnog-mapper on called ORFs, first optionally downloading the required databases
//

include { EGGNOG_DOWNLOAD } from '../../../modules/local/eggnog/download/main'
include { EGGNOG_MAPPER   } from '../../../modules/local/eggnog/mapper/main'
include { EGGNOG_SUM      } from '../../../modules/local/eggnog/sum/main'

workflow EGGNOG {
    take:
    faa
    collect_fcs

    main:

    EGGNOG_DOWNLOAD()

    ch_eggnog_database = EGGNOG_DOWNLOAD.out.eggnog_db
        .combine(EGGNOG_DOWNLOAD.out.dmnd)
        .combine(EGGNOG_DOWNLOAD.out.taxa_db)
        .combine(EGGNOG_DOWNLOAD.out.pkl)

    EGGNOG_MAPPER(faa, ch_eggnog_database)

    EGGNOG_SUM(EGGNOG_MAPPER.out.emappertsv, collect_fcs)

    emit:
    hits       = EGGNOG_MAPPER.out.hits
    emappertsv = EGGNOG_MAPPER.out.emappertsv
    sumtable   = EGGNOG_SUM.out.eggnog_summary
}
