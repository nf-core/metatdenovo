//
// Call ORFs with MetaEuk on contigs split by size, then concatenate output and reformat
// its GFF for downstream featureCounts
//

include { METAEUK_DOWNLOAD                    } from '../../../modules/local/metaeuk/download/main'
include { METAEUK_EASYPREDICT                 } from '../../../modules/nf-core/metaeuk/easypredict/main'
include { FIND_CONCATENATE as METAEUK_FAA_CAT } from '../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as METAEUK_GFF_CAT } from '../../../modules/nf-core/find/concatenate/main'
include { FORMAT_METAEUKFAA                   } from '../../../modules/local/format/metaeukfaa/main'
include { FORMAT_METAEUK_GFF                  } from '../../../modules/local/format/metaeuk/main'

workflow METAEUK {

    take:
    fasta     // channel: [ val(meta), path(fasta) ]
    db_path   // value: path to a pre-built database (file or directory), or falsy to auto-download
    db_name   // value: database name `metaeuk databases` understands, used only when db_path is falsy
    batchsize // channel: strings like '10.MB'. Usually from params.metaeuk_batchsize

    main:

    if ( db_path ) {
        ch_database = file(db_path, checkIfExists: true)
    } else {
        METAEUK_DOWNLOAD( db_name )
        ch_database = METAEUK_DOWNLOAD.out.database
    }

    // Calling is per-contig, so batching is safe and bounds extractorfs's peak memory.
    METAEUK_EASYPREDICT (
        fasta
            .map { _meta, contigs -> contigs }
            .splitFasta(size: batchsize, file: true)
            .map { ctg -> [ [ id: ctg.getBaseName() ], ctg ] },
        ch_database
    )

    ch_faa = fasta.map { meta, _contigs -> meta }
        .combine(METAEUK_EASYPREDICT.out.faa.collect { _meta, faa -> faa }.map { faas -> [ faas ] })

    METAEUK_FAA_CAT(ch_faa)

    ch_gff = fasta.map { meta, _contigs -> meta }
        .combine(METAEUK_EASYPREDICT.out.gff.collect { _meta, gff -> gff }.map { gffs -> [ gffs ] })

    METAEUK_GFF_CAT(ch_gff)

    // MetaEuk's fasta ids and gff ID= attributes disagree; both are reformatted to match.
    FORMAT_METAEUKFAA ( METAEUK_FAA_CAT.out.file_out )
    FORMAT_METAEUK_GFF ( METAEUK_GFF_CAT.out.file_out )

    emit:
    faa = FORMAT_METAEUKFAA.out.format_faa // channel: [ val(meta), path(faa) ]
    gff = FORMAT_METAEUK_GFF.out.format_gff // channel: [ val(meta), path(gff) ]
}
