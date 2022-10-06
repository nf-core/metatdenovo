//
// Run PROKKA on contigs that are split by size, then concatenate output and gunzip it
//

include { PROKKA               } from '../../modules/nf-core/prokka/main'
include { CAT_CAT as GFF_CAT   } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FAA_CAT   } from '../../modules/nf-core/cat/cat/main'
include { CAT_CAT as FNA_CAT   } from '../../modules/nf-core/cat/cat/main'

workflow PROKKA_CAT {
    take:
        contigs

    main:
        ch_versions = Channel.empty()

        PROKKA  (contigs.splitFasta(size: 10.MB, file: true).map { contigs -> [[id: contigs.getBaseName()], contigs] }, [], [] )
        ch_versions = ch_versions.mix(PROKKA.out.versions)
        GFF_CAT (PROKKA.out.gff.collect().map { [ [ id: 'prokka.gff.gz' ], it[1] ] } )
        FAA_CAT (PROKKA.out.faa.collect().map { [ [ id: 'prokka.faa.gz' ], it[1] ] } )
        FNA_CAT (PROKKA.out.fna.collect().map { [ [ id: 'prokka.fna.gz' ], it[1] ] } )
        ch_versions = ch_versions.mix(FNA_CAT.out.versions)

    emit:
        gff      = GFF_CAT.out.file_out
        faa      = FAA_CAT.out.file_out
        fna      = FNA_CAT.out.file_out

        versions = ch_versions

}
