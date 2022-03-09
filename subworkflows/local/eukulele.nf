//
// Run EUKULELE on contigs that are moved into a directory
//

include { EUKULELE        } from '../../modules/local/eukulele/main'
include { MV_DIR } from '../../modules/local/mvdir.nf'

workflow SUB_EUKULELE {
    take:
        contigs
        db


    main:
        ch_versions = Channel.empty()
        ch_eukulele_pathdb = Channel.empty()
        MV_DIR          ( contigs )
        EUKULELE        ( MV_DIR.out.contigs_dir, db )
        ch_versions = ch_versions.mix(EUKULELE.out.versions)
    
    emit:
        taxonomy_extimation = EUKULELE.out.taxonomy_extimation
        taxonomy_counts     = EUKULELE.out.taxonomy_counts
        diamond             = EUKULELE.out.diamond
        
        versions            = ch_versions
}
