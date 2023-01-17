//
// Run PROKKA on contigs that are split by size, then concatenate output and gunzip it
//

include { PROKKA               } from '../../modules/nf-core/prokka/main'
include { CAT_CAT as GFF_CAT   } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FAA_CAT   } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FNA_CAT   } from '../../modules/nf-core/cat/cat/main'

workflow PROKKA_SUBSETS {
    take:
        contigs // channel: [ val(meta), file(contigs) ]

    main:
        ch_versions = Channel.empty()

        contigs
            .splitFasta(size: 1.MB, file: true)
            .set { ch_prokka }
        PROKKA ( ch_prokka, [], [] )
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        GFF_CAT ( PROKKA.out.gff.collect().map { [ [ id: 'prokka' ], it[1] ] } )
        ch_versions = ch_versions.mix(GFF_CAT.out.versions)
        FAA_CAT ( PROKKA.out.faa.collect().map { [ [ id: 'prokka' ], it[1] ] } )
        ch_versions = ch_versions.mix(FAA_CAT.out.versions)
        FNA_CAT ( PROKKA.out.fna.collect().map { [ [ id: 'prokka' ], it[1] ] } )
        ch_versions = ch_versions.mix(FNA_CAT.out.versions)

    emit:
        gff      = GFF_CAT.out.file_out
        faa      = FAA_CAT.out.file_out
        fna      = FNA_CAT.out.file_out

        versions = ch_versions

}
