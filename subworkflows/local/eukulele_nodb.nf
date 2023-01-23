//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE                          } from '../../modules/local/eukulele/main'
include { FORMAT_TAX                        } from '../../modules/local/format_tax'

workflow SUB_EUKULELE_NODB {

    take:
        eukulele // Channel: val(meta), path(fasta), val(database), path(directory) 

    main:
        ch_versions = Channel.empty()
        EUKULELE( eukulele )
        FORMAT_TAX( EUKULELE.out.taxonomy_estimation.map { [ it[2], it[1] ] } )

    emit:
        taxonomy_estimation = EUKULELE.out.taxonomy_estimation
        taxonomy_counts     = EUKULELE.out.taxonomy_counts
        diamond             = EUKULELE.out.diamond
        versions            = EUKULELE.out.versions
}
