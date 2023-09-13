//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE_SEARCH                   } from '../../modules/local/eukulele/search'
include { EUKULELE_DOWNLOAD                 } from '../../modules/local/eukulele/download'
include { FORMAT_TAX                        } from '../../modules/local/format_tax'
include { SUM_TAXONOMY                      } from '../../modules/local/sum_taxonomy'

workflow SUB_EUKULELE {

    take:
        eukulele // Channel: val(meta), path(fasta), val(database), path(directory) 
        feature_counts

    main:
        ch_versions = Channel.empty()


        EUKULELE_DOWNLOAD ( eukulele.filter { it[2] }.map { [ it[2], it[3] ] } )
        ch_download = EUKULELE_DOWNLOAD.out.db

        Channel.empty()
            .mix ( EUKULELE_DOWNLOAD.out.db )
            .mix(eukulele.filter{ ! it[2] }.map { [ [], it[3] ] } )
            .merge( eukulele.map{ [ it[0], it[1] ] } )
            .map { [ [ id: "${it[2].id}.${it[0]}" ], it[3], it[0], it[1] ] } 
            .set { ch_eukulele }
        EUKULELE_SEARCH( ch_eukulele )

        FORMAT_TAX( EUKULELE_SEARCH.out.taxonomy_estimation.map { [ it[0], it[1] ] } )

        FORMAT_TAX.out.tax
            .join(ch_eukulele)
            .map { [ it[0], it[3], it[1] ] }
            .set { ch_sum_taxonomy }

        SUM_TAXONOMY ( ch_sum_taxonomy, feature_counts )

    emit:
        taxonomy_summary    = SUM_TAXONOMY.out.taxonomy_summary
        taxonomy_estimation = EUKULELE_SEARCH.out.taxonomy_estimation
        taxonomy_counts     = EUKULELE_SEARCH.out.taxonomy_counts
        diamond             = EUKULELE_SEARCH.out.diamond
        versions            = EUKULELE_SEARCH.out.versions
}
