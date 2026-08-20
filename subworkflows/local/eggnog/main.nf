//
// Run eggnog-mapper on called ORFs, first optionally downloading the required databases
//

include { EGGNOG_DOWNLOAD } from '../../../modules/local/eggnog/download/main'
include { EGGNOGMAPPER    } from '../../../modules/nf-core/eggnogmapper/main'
include { EGGNOG_FORMAT   } from '../../../modules/local/eggnog/format/main'
include { EGGNOG_SUM      } from '../../../modules/local/eggnog/sum/main'

workflow EGGNOG {
    take:
    faa           // channel: [ val(meta), path(faa) ]
    feature_counts // channel: [ val(meta), path(fcs) ] -- meta.caller must match faa's

    main:

    EGGNOG_DOWNLOAD()

    ch_search_mode_db = EGGNOG_DOWNLOAD.out.dmnd.map { dmnd -> [ 'diamond', dmnd ] }

    // The official module wants a single "eggnog_data" directory (--data_dir); stage the
    // db/taxa/pkl files into one at task-staging time only, via the module's own stageAs,
    // rather than restructuring EGGNOG_DOWNLOAD's own (storeDir-cached, so expensive to
    // invalidate) flat output layout.
    ch_eggnog_data_dir = EGGNOG_DOWNLOAD.out.eggnog_db
        .combine(EGGNOG_DOWNLOAD.out.taxa_db)
        .combine(EGGNOG_DOWNLOAD.out.pkl)

    EGGNOGMAPPER(faa, ch_search_mode_db, ch_eggnog_data_dir)

    EGGNOG_FORMAT(EGGNOGMAPPER.out.annotations)

    ch_eggnog_sum_input = EGGNOG_FORMAT.out.emappertsv
        .map { meta, tsv -> [ meta.caller, meta, tsv ] }
        .join( feature_counts.map { meta, fcs -> [ meta.caller, fcs ] } )
        .map { _caller, meta, tsv, fcs -> [ meta, tsv, fcs ] }

    EGGNOG_SUM(ch_eggnog_sum_input)

    emit:
    hits       = EGGNOGMAPPER.out.hits
    emappertsv = EGGNOG_FORMAT.out.emappertsv
    sumtable   = EGGNOG_SUM.out.eggnog_summary
}
