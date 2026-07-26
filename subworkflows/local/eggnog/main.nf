//
// Run eggnog-mapper on called ORFs, first optionally downloading the required databases
//

include { EGGNOG_DOWNLOAD }              from '../../../modules/local/eggnog/download/main'
include { EGGNOGMAPPER as EGGNOG_MAPPER } from '../../../modules/nf-core/eggnogmapper/main'
include { EGGNOG_FORMAT }                from '../../../modules/local/eggnog/format/main'
include { EGGNOG_SUM }                   from '../../../modules/local/eggnog/sum/main'

workflow EGGNOG {
    take:
    faa
    collect_fcs

    main:

    EGGNOG_DOWNLOAD()

    ch_search_mode_db = EGGNOG_DOWNLOAD.out.dmnd.map { dmnd -> [ 'diamond', dmnd ] }

    EGGNOG_MAPPER(faa, ch_search_mode_db, EGGNOG_DOWNLOAD.out.eggnog_data_dir)

    EGGNOG_FORMAT(EGGNOG_MAPPER.out.annotations)

    EGGNOG_SUM(EGGNOG_FORMAT.out.emappertsv, collect_fcs)

    emit:
    hits       = EGGNOG_MAPPER.out.hits
    emappertsv = EGGNOG_FORMAT.out.emappertsv
    sumtable   = EGGNOG_SUM.out.eggnog_summary
}
