//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUK_SRUN                         } from '../../modules/local/eukulele/main_srun'
include { EUKULELE_SDB                     } from '../../modules/local/eukulele/download_second_db'
include { FORMAT_TAX as FORMAT_SECOND_TAX  } from '../../modules/local/format_tax'

workflow EUKULELE_SRUN {

    take:
        fastaprot

    main:
        ch_versions = Channel.empty()

        String directoryName = params.eukulele_sdbpath
        File directory = new File(directoryName)
        String eukuleleDB = params.eukulele_second_db
        File test = new File(eukuleleDB)

        if ( params.eukulele_multirun ) {
            if(! directory.exists() ) {
                directory.mkdir()
                EUKULELE_SDB( )
                if ( params.eukulele_second_db == 'mmetsp' ) {
                    EUK_SRUN(fastaprot, EUKULELE_SDB.out.mmetsp_db)
                    FORMAT_SECOND_TAX(EUK_SRUN.out.taxonomy_estimation.map { it[1] } )
                } else
                    EUK_SRUN(fastaprot, EUKULELE_SDB.out.phylo_db)
                    FORMAT_SECOND_TAX(EUK_SRUN.out.taxonomy_estimation.map { it[1] } )
            } else {
                    ch_database = Channel.fromPath(params.eukulele_sdbpath)
                    EUK_SRUN(fastaprot, ch_database)
                    FORMAT_SECOND_TAX(EUK_SRUN.out.taxonomy_estimation.map { it[1] } )
            }
        }

    emit:
        taxonomy_estimation = EUK_SRUN.out.taxonomy_estimation
        taxonomy_counts     = EUK_SRUN.out.taxonomy_counts
        diamond             = EUK_SRUN.out.diamond
        versions            = EUK_SRUN.out.versions
}
