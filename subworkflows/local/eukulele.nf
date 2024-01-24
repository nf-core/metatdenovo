//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE_SEARCH                   } from '../../modules/local/eukulele/search'
include { FORMAT_TAX                        } from '../../modules/local/format_tax'
include { SUM_TAXONOMY                      } from '../../modules/local/sum_taxonomy'

workflow SUB_EUKULELE {

    take:
    eukulele // Channel: val(meta), path(fasta), val(database), path(directory)
    feature_counts

    main:
    ch_versions = Channel.empty()

    EUKULELE_SEARCH( eukulele )
    ch_versions = ch_versions.mix ( EUKULELE_SEARCH.out.versions )

    FORMAT_TAX( EUKULELE_SEARCH.out.taxonomy_estimation.map { meta, taxonomy, dbname -> [ meta, taxonomy ] } )
    ch_versions = ch_versions.mix ( FORMAT_TAX.out.versions )

    FORMAT_TAX.out.tax
        .join(eukulele)
        .map { meta, taxonomy, protein, dbname, database -> [ meta, dbname, taxonomy ] }
        .set { ch_sum_taxonomy }

    SUM_TAXONOMY ( ch_sum_taxonomy, feature_counts )
    ch_versions = ch_versions.mix ( SUM_TAXONOMY.out.versions )

    emit:
    taxonomy_summary    = SUM_TAXONOMY.out.taxonomy_summary
    taxonomy_estimation = EUKULELE_SEARCH.out.taxonomy_estimation
    taxonomy_counts     = EUKULELE_SEARCH.out.taxonomy_counts
    diamond             = EUKULELE_SEARCH.out.diamond
    versions            = ch_versions
}
