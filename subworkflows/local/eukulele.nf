//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE                          } from '../../modules/local/eukulele/main'
include { EUKULELE_DOWNLOAD                 } from '../../modules/local/eukulele/download'
include { FORMAT_TAX                        } from '../../modules/local/format_tax'

workflow SUB_EUKULELE {

    take:
        eukulele // Channel: val(meta), path(fasta), val(database), path(directory) 

    main:
        ch_versions = Channel.empty()
        EUKULELE_DOWNLOAD ( eukulele.map { it[3] }, eukulele.map { it[2] } )
        eukulele
            .map { [ it[0], it[1], it[2] ] }
            .combine( EUKULELE_DOWNLOAD.out.db )
            .set { ch_eukulele }
        EUKULELE( ch_eukulele )

        FORMAT_TAX( EUKULELE.out.taxonomy_estimation.map { [ it[2], it[1] ] } )

    emit:
        taxonomy_estimation = EUKULELE.out.taxonomy_estimation
        taxonomy_counts     = EUKULELE.out.taxonomy_counts
        diamond             = EUKULELE.out.diamond
        versions            = EUKULELE.out.versions
}
