//
// Run PROKKA on contigs that are split by size, then concatenate output and gunzip it
//

include { PROKKA               } from '../../modules/nf-core/prokka/main'
include { CAT_CAT as GFF_CAT   } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FAA_CAT   } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FNA_CAT   } from '../../modules/nf-core/cat/cat/main'
include { PROKKAGFF2TSV        } from '../../modules/local/prokkagff2tsv'

workflow PROKKA_SUBSETS {
    take:
        contigs // channel:  tuple val(meta), file(contigs)

    main:
        ch_versions = Channel.empty()

        PROKKA ( contigs.map{ it[1] }.splitFasta(size: 10.MB, file: true).map { contigs -> [ [ id: contigs.getBaseName() ], contigs] }, [], []  )
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        ch_log  = PROKKA.out.txt.collect()

        contigs.map{ [ id:"${it[0].id}.prokka" ] }
            .combine(PROKKA.out.gff.collect { it[1] }.map { [ it ] })
            .set { ch_gff }
        GFF_CAT ( ch_gff )
        ch_versions = ch_versions.mix(GFF_CAT.out.versions)

        contigs.map{ [ id:"${it[0].id}.prokka" ] }
            .combine(PROKKA.out.faa.collect { it[1] }.map { [ it ] })
            .set { ch_faa }
        FAA_CAT ( ch_faa )
        ch_versions = ch_versions.mix(FAA_CAT.out.versions)

        contigs.map{ [ id:"${it[0].id}.prokka" ] }
            .combine(PROKKA.out.fna.collect { it[1] }.map { [ it ] })
            .set { ch_fna }
        FNA_CAT ( ch_fna )
        ch_versions = ch_versions.mix(FNA_CAT.out.versions)

        PROKKAGFF2TSV (
            GFF_CAT.out.file_out
        )
        ch_versions = ch_versions.mix(PROKKAGFF2TSV.out.versions)

    emit:
        gff        = GFF_CAT.out.file_out
        faa        = FAA_CAT.out.file_out
        fna        = FNA_CAT.out.file_out
        //gfftsv     = PROKKAGFF2TSV.out.tsv
        prokka_log = ch_log
        versions   = ch_versions

}
