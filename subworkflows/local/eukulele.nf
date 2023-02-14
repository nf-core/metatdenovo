//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE                          } from '../../modules/local/eukulele/main'
include { EUKULELE_DOWNLOAD                 } from '../../modules/local/eukulele/download'
include { FORMAT_TAX                        } from '../../modules/local/format_tax'
include { SUM_TAXONOMY                      } from '../../modules/local/sum_taxonomy'

workflow SUB_EUKULELE {

    take:
        eukulele // Channel: val(meta), path(fasta), val(database), path(directory) 
        collect_fcs

    main:
        ch_versions = Channel.empty()
        EUKULELE_DOWNLOAD ( eukulele.filter{ it[2] }.map { [ it[2], it[3] ] } )
        ch_download = EUKULELE_DOWNLOAD.out.db
        Channel.empty()
            .mix ( EUKULELE_DOWNLOAD.out.db )
            .mix(eukulele.filter{ ! it[2] }.map { [ [], it[3] ] } )
            .merge( eukulele.map{ [ it[0], it[1] ] } )
            .map { [ it[2], it[3], it[0], it[1] ] } 
            .set { ch_eukulele }
        EUKULELE( ch_eukulele )

        FORMAT_TAX( EUKULELE.out.taxonomy_estimation.map { [ [id: it[2]], it[1] ] } )
        SUM_TAXONOMY( FORMAT_TAX.out.tax, collect_fcs )

    emit:
        taxonomy_summary    = SUM_TAXONOMY.out.taxonomy_summary
        taxonomy_estimation = EUKULELE.out.taxonomy_estimation
        taxonomy_counts     = EUKULELE.out.taxonomy_counts
        diamond             = EUKULELE.out.diamond
        versions            = EUKULELE.out.versions
}
