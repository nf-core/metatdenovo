//
// Run eggnog-mapper on called ORFs, first optionally downloading the required databases
//

include { EGGNOG_DOWNLOAD } from '../../modules/local/eggnog/download'
include { EGGNOG_MAPPER   } from '../../modules/local/eggnog/mapper'
include { EGGNOG_TABLE    } from '../../modules/local/eggnog_table'

workflow EGGNOG {
    take:
        faa

    main:
        ch_versions = Channel.empty()

        String directoryName = params.eggnog_dbpath
        File directory = new File(directoryName)
        String eggnogDB = params.eggnog_dbpath + "eggnog.db"
        File test = new File(eggnogDB)

        if (! directory.exists()){
            EGGNOG_DOWNLOAD()
            EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
            EGGNOG_TABLE(EGGNOG_MAPPER.out.annotations)
        } else {
            if (! test.exists()){
                EGGNOG_DOWNLOAD()
                EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
                EGGNOG_TABLE(EGGNOG_MAPPER.out.annotations)
            }  else {
            ch_dbpath = Channel.fromPath(params.eggnog_dbpath)
            EGGNOG_MAPPER(faa, ch_dbpath)
            EGGNOG_TABLE(EGGNOG_MAPPER.out.annotations)
            }
        }

        ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)
        ch_versions = ch_versions.mix(EGGNOG_TABLE.out.versions)

    emit:
        hits     = EGGNOG_MAPPER.out.hits
        versions = ch_versions

}
