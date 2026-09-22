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
    db_url         // string: URL of eggnog.db.gz
    dmnd_url       // string: URL of eggnog_proteins.dmnd.gz
    taxa_url       // string: URL of eggnog.taxa.tar.gz

    main:

    EGGNOG_DOWNLOAD( file(db_url), file(dmnd_url), file(taxa_url) )

    ch_search_mode_db = EGGNOG_DOWNLOAD.out.dmnd.map { dmnd -> [ 'diamond', dmnd ] }

    // EGGNOGMAPPER wants one data dir; stageAs builds it, leaving EGGNOG_DOWNLOAD's storeDir layout alone.
    ch_eggnog_data_dir = EGGNOG_DOWNLOAD.out.eggnog_db
        .combine(EGGNOG_DOWNLOAD.out.taxa_db)
        .combine(EGGNOG_DOWNLOAD.out.pkl)
        // .combine() demotes to a one-item queue channel, which would cap EGGNOGMAPPER at one task.
        .first()

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
