//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN     } from '../../modules/local/eukulele/kofamscan'
include { DOWNLOAD_KODB } from '../../modules/local/eukulele/download_kofamscan_db'

workflow SUB_KOFAMSCAN {

    take:
        kofamscan // Channel: val(meta), path(fasta), val(database), path(ko_list), path(koprofiles)

    main:
        ch_versions = Channel.empty()
        DOWNLOAD_KODB ( eukulele.filter{ it[2] }.map { [ it[2], it[3] ] } )
        ch_download = EUKULELE_DOWNLOAD.out.db
        Channel.empty()
            .mix ( EUKULELE_DOWNLOAD.out.db )
            .mix(eukulele.filter{ ! it[2] }.map { [ [], it[3] ] } )
            .merge( eukulele.map{ [ it[0], it[1] ] } )
            .map { [ it[2], it[3], it[0], it[1] ] } 
            .set { ch_eukulele }
        EUKULELE( ch_eukulele )

    emit:
        taxonomy_summary    = SUM_TAXONOMY.out.taxonomy_summary
        versions            = EUKULELE.out.versions
}
