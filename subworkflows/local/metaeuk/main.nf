//
// Call ORFs with MetaEuk on contigs split by size, then concatenate output and reformat
// its GFF for downstream featureCounts
//

include { METAEUK_EASYPREDICT                 } from '../../../modules/nf-core/metaeuk/easypredict/main'
include { FIND_CONCATENATE as METAEUK_FAA_CAT } from '../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as METAEUK_GFF_CAT } from '../../../modules/nf-core/find/concatenate/main'
include { FORMAT_METAEUKFAA                   } from '../../../modules/local/format/metaeukfaa/main'
include { FORMAT_METAEUK_GFF                  } from '../../../modules/local/format/metaeuk/main'

workflow METAEUK {

    take:
    fasta     // channel: [ val(meta), path(fasta) ]
    database  // channel: path(database)
    batchsize // channel: strings like '10.MB'. Usually from params.metaeuk_batchsize

    main:

    // MetaEuk's splice-aware ORF calling is per-contig, so splitting the assembly into
    // batches (same rationale/mechanism as PROKKA_SUBSETS) is safe, and directly bounds
    // extractorfs's peak memory, which scales with however many contigs it's handed at
    // once -- unlike Prokka, a killed/aborted batch here only costs re-running that batch,
    // not the whole multi-hour, whole-assembly task.
    METAEUK_EASYPREDICT (
        fasta
            .map { _meta, contigs -> contigs }
            .splitFasta(size: batchsize, file: true)
            .map { ctg -> [ [ id: ctg.getBaseName() ], ctg ] },
        database
    )

    ch_faa = fasta.map { meta, _contigs -> meta }
        .combine(METAEUK_EASYPREDICT.out.faa.collect { _meta, faa -> faa }.map { faas -> [ faas ] })

    METAEUK_FAA_CAT(ch_faa)

    ch_gff = fasta.map { meta, _contigs -> meta }
        .combine(METAEUK_EASYPREDICT.out.gff.collect { _meta, gff -> gff }.map { gffs -> [ gffs ] })

    METAEUK_GFF_CAT(ch_gff)

    // Both outputs get reformatted so the protein fasta's sequence ids and the gff's ID= attributes
    // are the same string -- MetaEuk's own two outputs disagree, see FORMAT_METAEUKFAA.
    FORMAT_METAEUKFAA ( METAEUK_FAA_CAT.out.file_out )
    FORMAT_METAEUK_GFF ( METAEUK_GFF_CAT.out.file_out )

    emit:
    faa = FORMAT_METAEUKFAA.out.format_faa // channel: [ val(meta), path(faa) ]
    gff = FORMAT_METAEUK_GFF.out.format_gff // channel: [ val(meta), path(gff) ]
}
