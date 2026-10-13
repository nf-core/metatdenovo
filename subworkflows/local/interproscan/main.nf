//
// Run InterProScan on called ORFs, first downloading the required database
//

include { INTERPROSCAN_DATABASE } from '../../../modules/local/interproscan/database/main'
include { INTERPROSCAN          } from '../../../modules/nf-core/interproscan/main'
include { INTERPROSCAN_FORMAT   } from '../../../modules/local/interproscan/format/main'
include { INTERPROSCAN_SUM      } from '../../../modules/local/interproscan/sum/main'

workflow INTERPROSCAN_ANNOTATION {
    take:
    faa            // channel: [ val(meta), path(faa) ]
    feature_counts // channel: [ val(meta), path(fcs) ] -- meta.caller must match faa's
    db_url         // string: URL of the InterProScan data release tarball
    batchsize      // integer: residues per batch INTERPROSCAN is run on

    main:

    INTERPROSCAN_DATABASE( file(db_url) )

    // splitFasta cannot read a zero-byte gzip, which stub runs produce
    ch_interproscan_split = faa
        .branch { _meta, f ->
            empty: f.size() == 0
            other: true
        }

    INTERPROSCAN(
        ch_interproscan_split.other
            .splitFasta(size: batchsize, file: true, elem: 1)
            .map { meta, f -> [ meta, f, f.baseName.tokenize('.').last() as Integer ] }
            .mix( ch_interproscan_split.empty.map { meta, f -> [ meta, f, 1 ] } )
            .map { meta, f, n -> [ meta + [ id: String.format('%s.%03d', meta.id, n), parent: meta ], f ] },
        INTERPROSCAN_DATABASE.out.db
    )

    // Sorted by batch number: input order is part of the downstream task hashes and the merged row order
    ch_interproscan_tsv = INTERPROSCAN.out.tsv
        .map { meta, tsv -> [ meta.parent, tsv ] }
        .groupTuple(sort: { a, b -> a.name <=> b.name })

    INTERPROSCAN_FORMAT(ch_interproscan_tsv)

    ch_interproscan_sum_input = INTERPROSCAN_FORMAT.out.tsv
        .map { meta, tsv -> [ meta.caller, meta, tsv ] }
        .join( feature_counts.map { meta, fcs -> [ meta.caller, fcs ] } )
        .map { _caller, meta, tsv, fcs -> [ meta, tsv, fcs ] }

    INTERPROSCAN_SUM(ch_interproscan_sum_input)

    emit:
    annotation = INTERPROSCAN_FORMAT.out.tsv
    sumtable   = INTERPROSCAN_SUM.out.interproscan_summary
}
