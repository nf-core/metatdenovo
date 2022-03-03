//
// Run PROKKA on contigs that are split by size, then concatenate output and gunzip it
//

include { EUKULELE        } from '../../modules/local/eukulele/main'
include { STAGE_FASTA_DIR } from '../../modules/local/mvdir.nf'

workflow SUB_EUKULELE {
    take:
        contigs

    main:
        ch_versions = Channel.empty()
        
        STAGE_FASTA_DIR ( contigs )
        EUKULELE        ( STAGE_FASTA_DIR.out.contigs_dir )
        ch_versions = ch_versions.mix(EUKULELE.out.versions)
    
    emit:
        taxonomy_extimation = EUKULELE.out.taxonomy_extimation
        taxonomy_counts     = EUKULELE.out.taxonomy_counts
        diamond             = EUKULELE.out.diamond
        
        versions            = ch_versions
}
