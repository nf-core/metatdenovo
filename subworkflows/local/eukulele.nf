//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE                          } from '../../modules/local/eukulele/main'
include { EUKULELE_DB                       } from '../../modules/local/eukulele/download'
include { FORMAT_TAX                        } from '../../modules/local/format_tax'

workflow SUB_EUKULELE {

    take:
        fastaprot
        eukulele_db

    main:
        ch_versions = Channel.empty()

        String directoryName = params.eukulele_dbpath
        File directory = new File(directoryName)
        if(! directory.exists() ) {
            directory.mkdir()
            EUKULELE_DB(eukulele_db)
            EUKULELE(fastaprot,EUKULELE_DB.out.db, eukulele_db)
            } else {
                ch_dbpath = Channel.fromPath(params.eukulele_dbpath)
                fastaprot
                    .combine(eukulele_db)
                    .combine(ch_dbpath)
                    .set { ch_eukulele }
                EUKULELE( ch_eukulele )
            }

            FORMAT_TAX( EUKULELE.out.taxonomy_estimation )


    emit:
        taxonomy_estimation = EUKULELE.out.taxonomy_estimation
        taxonomy_counts     = EUKULELE.out.taxonomy_counts
        diamond             = EUKULELE.out.diamond
        versions            = EUKULELE.out.versions
}
