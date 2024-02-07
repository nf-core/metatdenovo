//
// Run eggnog-mapper on called ORFs, first optionally downloading the required databases
//

include { EGGNOG_DOWNLOAD } from '../../modules/local/eggnog/download'
include { EGGNOG_MAPPER   } from '../../modules/local/eggnog/mapper'
include { EGGNOG_SUM      } from '../../modules/local/eggnog/sum'

workflow EGGNOG {
    take:
    eggnog_dbpath
    faa
    collect_fcs

    main:
    ch_versions = Channel.empty()
    String eggnogDB      = eggnog_dbpath + "/eggnog.db"
    File eggnogfile      = new File(eggnogDB)

    if ( ! eggnogfile.exists() ) {
        EGGNOG_DOWNLOAD()
        ch_dbpath = EGGNOG_DOWNLOAD.out.db
        ch_versions = ch_versions.mix ( EGGNOG_DOWNLOAD.out.versions )
    } else {
        ch_dbpath = Channel.fromPath(eggnog_dbpath, checkIfExists: true)
    }

    EGGNOG_MAPPER ( faa, ch_dbpath)
    ch_versions = ch_versions.mix ( EGGNOG_MAPPER.out.versions )

    EGGNOG_SUM ( EGGNOG_MAPPER.out.emappertsv, collect_fcs )
    ch_versions = ch_versions.mix ( EGGNOG_SUM.out.versions )

    emit:
    hits       = EGGNOG_MAPPER.out.hits
    emappertsv = EGGNOG_MAPPER.out.emappertsv
    versions   = ch_versions
    sumtable   = EGGNOG_SUM.out.eggnog_summary

}
