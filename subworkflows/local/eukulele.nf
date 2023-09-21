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
        collect_fcs

    main:
        ch_versions = Channel.empty()

        String directoryName = params.eukulele_dbpath
        File directory       = new File(directoryName)
        // get files in euk directory, and checks if there is a reference.pep.fa in the
        // first one. Not the most robust but if this fails it will simply download
        // the database anyways.
        List euk_files       = []
        new File(directoryName).eachFile() {
            file-> euk_files.add(file)
        }
        String eukdb         = euk_files.get(0).toString()
        File eukpepfa        = new File(eukdb.plus("/reference.pep.fa"))

        if ( ! directory.exists() ) {
            directory.mkdir()
        }

        if ( ! eukpepfa.exists() ) {
            EUKULELE_DOWNLOAD ( eukulele.filter{ it[2] }.map { [ it[2], it[3] ] } )
            ch_download = EUKULELE_DOWNLOAD.out.db
        } else {
            // tuple val("${db}"), path("${directory}/${db}"), emit: db
            // mimic download output with subdir and db name
            String db = eukdb.split("/")[0]
            Channel.empty()
                .map { tuple val(db), Channel.fromPath(eukdb) }
                .set { ch_download }
        }

        Channel.empty()
            .mix ( ch_download )
            .mix(eukulele.filter{ ! it[2] }.map { [ [], it[3] ] } )
            .merge( eukulele.map{ [ it[0], it[1] ] } )
            .map { [ [ id: "${it[2].id}.${it[0]}" ], it[3], it[0], it[1] ] }
            .set { ch_eukulele }
        EUKULELE_SEARCH( ch_eukulele )

        FORMAT_TAX( EUKULELE_SEARCH.out.taxonomy_estimation.map { [ it[0], it[1] ] } )
        SUM_TAXONOMY( FORMAT_TAX.out.tax, collect_fcs )

    emit:
        taxonomy_summary    = SUM_TAXONOMY.out.taxonomy_summary
        taxonomy_estimation = EUKULELE_SEARCH.out.taxonomy_estimation
        taxonomy_counts     = EUKULELE_SEARCH.out.taxonomy_counts
        diamond             = EUKULELE_SEARCH.out.diamond
        versions            = EUKULELE_SEARCH.out.versions
}
