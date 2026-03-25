//
// Run PROKKA on contigs that are split by size, then concatenate output and gunzip it
//

include { PROKKA               } from '../../../../modules/nf-core/prokka/main'
include { CAT_CAT as GFF_CAT   } from '../../../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FAA_CAT   } from '../../../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FFN_CAT   } from '../../../../modules/nf-core/cat/cat/main'
include { PROKKAGFF2TSV        } from '../../../../modules/local/prokka/gff2tsv/main'

workflow PROKKA_SUBSETS {
    take:
    contigs   // channel:  tuple val(meta), file(contigs)
    batchsize // channel: strings like '10.MB'. Usually from params.prokka_batchsize

    main:
    ch_versions = channel.empty()

    PROKKA ( contigs.map{ _meta, ctg -> ctg }.splitFasta(size: batchsize, file: true).map { ctg -> [ [ id: ctg.getBaseName() ], ctg] }, [], []  )
    ch_versions = ch_versions.mix(PROKKA.out.versions)
    ch_log  = PROKKA.out.txt.map { _meta, log -> log }.collect()
    contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.gff.collect { _meta, gff -> gff }.map { gff -> [ gff ] })
        .set { ch_gff }
    GFF_CAT ( ch_gff )

    contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.faa.collect { _meta, protein -> protein }.map { protein -> [ protein ] })
        .set { ch_faa }
    FAA_CAT ( ch_faa )
    //ch_versions = ch_versions.mix(FAA_CAT.out.versions)

    contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.ffn.collect { _meta, fnn -> fnn }.map { fnn -> [ fnn ] })
        .set { ch_ffn }
    FFN_CAT ( ch_ffn )
    //ch_versions = ch_versions.mix(FFN_CAT.out.versions)

    PROKKAGFF2TSV ( GFF_CAT.out.file_out)
    //ch_versions = ch_versions.mix(PROKKAGFF2TSV.out.versions)

    emit:
    gff        = GFF_CAT.out.file_out.first()
    faa        = FAA_CAT.out.file_out.first()
    ffn        = FFN_CAT.out.file_out.first()
    gfftsv     = PROKKAGFF2TSV.out.tsv
    prokka_log = ch_log
    //versions   = ch_versions

}
