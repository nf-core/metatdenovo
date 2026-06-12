//
// Run PROKKA on contigs that are split by size, then concatenate output and gunzip it
//

include { PROKKA                      } from '../../../../modules/nf-core/prokka/main'
include { FIND_CONCATENATE as GFF_CAT } from '../../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as FAA_CAT } from '../../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as FFN_CAT } from '../../../../modules/nf-core/find/concatenate/main'
include { PROKKAGFF2TSV               } from '../../../../modules/local/prokka/gff2tsv/main'

workflow PROKKA_SUBSETS {
    take:
    contigs   // channel:  tuple val(meta), file(contigs)
    batchsize // channel: strings like '10.MB'. Usually from params.prokka_batchsize

    main:

    PROKKA(
        contigs
            .map { _meta, ctg -> ctg }
            .splitFasta(size: batchsize, file: true)
            .map { ctg -> [ [ id: ctg.getBaseName() ], ctg] },
        [],
        []
    )

    ch_log  = PROKKA.out.txt.map { _meta, log -> log }.collect()

    ch_gff = contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.gff.collect { _meta, gff -> gff }.map { gff -> [ gff ] })

    GFF_CAT(ch_gff)

    ch_faa = contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.faa.collect { _meta, protein -> protein }.map { protein -> [ protein ] })

    FAA_CAT(ch_faa)

    ch_ffn = contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.ffn.collect { _meta, fnn -> fnn }.map { fnn -> [ fnn ] })

    FFN_CAT(ch_ffn)

    PROKKAGFF2TSV(GFF_CAT.out.file_out)

    emit:
    gff        = GFF_CAT.out.file_out.first()
    faa        = FAA_CAT.out.file_out.first()
    ffn        = FFN_CAT.out.file_out.first()
    gfftsv     = PROKKAGFF2TSV.out.tsv
    prokka_log = ch_log
}
