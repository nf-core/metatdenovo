//
// Call ORFs with MetaEuk, then reformat its GFF for downstream featureCounts
//

include { METAEUK_EASYPREDICT } from '../../../modules/nf-core/metaeuk/easypredict/main'
include { FORMAT_METAEUK_FAA  } from '../../../modules/local/format/metaeuk_faa/main'
include { FORMAT_METAEUK_GFF  } from '../../../modules/local/format/metaeuk/main'

workflow METAEUK {

    take:
    fasta    // channel: [ val(meta), path(fasta) ]
    database // channel: path(database)

    main:

    METAEUK_EASYPREDICT ( fasta, database )

    // Both outputs get reformatted so the protein fasta's sequence ids and the gff's ID= attributes
    // are the same string -- MetaEuk's own two outputs disagree, see FORMAT_METAEUK_FAA.
    FORMAT_METAEUK_FAA ( METAEUK_EASYPREDICT.out.faa )
    FORMAT_METAEUK_GFF ( METAEUK_EASYPREDICT.out.gff )

    emit:
    faa = FORMAT_METAEUK_FAA.out.format_faa // channel: [ val(meta), path(faa) ]
    gff = FORMAT_METAEUK_GFF.out.format_gff // channel: [ val(meta), path(gff) ]
}
