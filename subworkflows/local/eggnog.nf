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

        EGGNOG_DOWNLOAD()
        EGGNOG_MAPPER(faa, EGGNOG_DOWNLOAD.out.db)

        ch_versions = ch_versions.mix(EGGNOG_MAPPER.out.versions)

    emit:
        hits     = EGGNOG_MAPPER.out.hits
        versions = ch_versions

}
