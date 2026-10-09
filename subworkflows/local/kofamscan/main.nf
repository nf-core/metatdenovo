//
// Run KOFAMSCAN on protein fasta from orf_caller output
//

include { KOFAMSCAN_DOWNLOAD }         from '../../../modules/local/kofamscan/download/main'
include { KOFAMSCAN as KOFAMSCAN_SCAN } from '../../../modules/nf-core/kofamscan/main'
include { KOFAMSCAN_FORMAT }           from '../../../modules/local/kofamscan/format/main'
include { KOFAMSCAN_UNIQUE }           from '../../../modules/local/kofamscan/unique/main'
include { KOFAMSCAN_SUM }              from '../../../modules/local/kofamscan/sum/main'

workflow KOFAMSCAN {

    take:
    kofamscan       // Channel: val(meta), path(fasta)
    fcs             // channel: [ val(meta), path(fcs) ] -- meta.caller must match kofamscan's
    ko_list_url     // string: URL to download the KOfam ko_list file from
    profiles_url    // string: URL to download the KOfam HMM profiles archive from
    batchsize       // integer: residues per KofamScan batch (splitFasta size counts sequence characters only)

    main:

    KOFAMSCAN_DOWNLOAD( file(ko_list_url), file(profiles_url) )

    // splitFasta cannot read a zero-byte gzip, which stub runs produce
    ch_kofamscan_split = kofamscan
        .branch { _meta, faa ->
            empty: faa.size() == 0
            other: true
        }

    // E-values are relative to each batch's size; KO assignment uses score thresholds only
    KOFAMSCAN_SCAN(
        ch_kofamscan_split.other
            .splitFasta(size: batchsize, file: true, elem: 1)
            .map { meta, faa -> [ meta, faa, faa.baseName.tokenize('.').last() as Integer ] }
            .mix( ch_kofamscan_split.empty.map { meta, faa -> [ meta, faa, 1 ] } )
            .map { meta, faa, n -> [ meta + [ id: String.format('%s.%03d', meta.id, n), parent: meta ], faa ] },
        KOFAMSCAN_DOWNLOAD.out.koprofiles,
        KOFAMSCAN_DOWNLOAD.out.ko_list
    )

    // sort: input order is part of the downstream task hashes
    ch_kofamscan_tsv = KOFAMSCAN_SCAN.out.tsv
        .map { meta, tsv -> [ meta.parent, tsv ] }
        .groupTuple(sort: true)

    KOFAMSCAN_FORMAT( ch_kofamscan_tsv )

    KOFAMSCAN_UNIQUE( KOFAMSCAN_FORMAT.out.kofamtsv )

    ch_kofamscan_sum_input = ch_kofamscan_tsv
        .map { meta, tsv -> [ meta.caller, meta, tsv ] }
        .join( fcs.map { meta, f -> [ meta.caller, f ] } )
        .map { _caller, meta, tsv, f -> [ meta, tsv, f ] }

    KOFAMSCAN_SUM( ch_kofamscan_sum_input )

    emit:
    kofam_table_out   = ch_kofamscan_tsv
    kofam_table_tsv   = KOFAMSCAN_FORMAT.out.kofamtsv
    kofam_table_uniq  = KOFAMSCAN_UNIQUE.out.kofamuniq
    kofamscan_summary = KOFAMSCAN_SUM.out.kofamscan_summary
}
