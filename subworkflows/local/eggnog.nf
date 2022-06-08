//
// Run eggnog-mapper on called ORFs, first optionally downloading the required databases
//

include { EGGNOG_DOWNLOAD } from '../../modules/local/eggnog/download'
include { EGGNOG_MAPPER   } from '../../modules/local/eggnog/mapper'

workflow EGGNOG {
    take:
        faa

    main:
        ch_versions = Channel.empty()

        String directoryName = params.eggnog_dbpath
        File directory = new File(directoryName)
        if (! directory.exists()){
            EGGNOG_DOWNLOAD()
            EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)
        } else {
            ch_dbpath = Channel.fromPath(params.eggnog_dbpath)
            EGGNOG_MAPPER(faa, ch_dbpath)
        }

        ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)

    emit:
        hits     = EGGNOG_MAPPER.out.hits
        versions = ch_versions

}
