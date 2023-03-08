//
// Run eggnog-mapper on called ORFs, first optionally downloading the required databases
//

include { EGGNOG_DOWNLOAD } from '../../modules/local/eggnog/download'
include { EGGNOG_MAPPER   } from '../../modules/local/eggnog/mapper'
include { EGGNOG_TABLE    } from '../../modules/local/eggnog_table'
include { SUM_EGGNOG      } from '../../modules/local/sum_eggnog.nf'

workflow EGGNOG {
    take:
        faa
        collect_fcs
    main:
        ch_versions = Channel.empty()

        String directoryName = params.eggnog_dbpath
        File directory = new File(directoryName)
        String eggnogDB = params.eggnog_dbpath + "eggnog.db"
        File test = new File(eggnogDB)

        if (! directory.exists()){
            directory.mkdir()
            EGGNOG_DOWNLOAD()
            EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
            EGGNOG_TABLE(EGGNOG_MAPPER.out.annotations)
            SUM_EGGNOG(EGGNOG_TABLE.out.eggtab, collect_fcs )
        } else {
            if (! test.exists()){
                EGGNOG_DOWNLOAD()
                EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
                EGGNOG_TABLE(EGGNOG_MAPPER.out.annotations)
                SUM_EGGNOG(EGGNOG_TABLE.out.eggtab, collect_fcs )
            }  else {
            ch_dbpath = Channel.fromPath(params.eggnog_dbpath)
            EGGNOG_MAPPER(faa, ch_dbpath)
            EGGNOG_TABLE(EGGNOG_MAPPER.out.annotations)
            SUM_EGGNOG(EGGNOG_TABLE.out.eggtab, collect_fcs )
            }
        }

        //ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)
        ch_versions = ch_versions.mix(EGGNOG_TABLE.out.versions)
        ch_versions = ch_versions.mix(SUM_EGGNOG.out.versions)

    emit:
        hits              = EGGNOG_MAPPER.out.hits
        formateggnogtable = EGGNOG_TABLE.out.eggtab
        versions          = ch_versions
        sumtable          = SUM_EGGNOG.out.eggnog_summary

}
