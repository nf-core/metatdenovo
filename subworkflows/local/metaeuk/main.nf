//
// Call ORFs with MetaEuk, then reformat its GFF for downstream featureCounts
//

include { METAEUK_EASYPREDICT } from '../../../modules/nf-core/metaeuk/easypredict/main'
include { FORMAT_METAEUK_GFF  } from '../../../modules/local/format/metaeuk/main'

workflow METAEUK {

    take:
    fasta    // channel: [ val(meta), path(fasta) ]
    database // channel: path(database)

    main:

    METAEUK_EASYPREDICT ( fasta, database )

    FORMAT_METAEUK_GFF ( METAEUK_EASYPREDICT.out.gff )

    emit:
    faa = METAEUK_EASYPREDICT.out.faa       // channel: [ val(meta), path(faa) ]
    gff = FORMAT_METAEUK_GFF.out.format_gff // channel: [ val(meta), path(gff) ]
}
