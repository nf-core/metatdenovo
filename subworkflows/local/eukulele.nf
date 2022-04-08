//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE      } from '../../modules/local/eukulele/main'
include { MV_DIR        } from '../../modules/local/mvdir'
include { EUKULELE_DB   } from '../../modules/local/eukulele/download'

workflow SUB_EUKULELE {

    take:
        fastaprot

    main:
        
        EUKULELE_DB( )
        MV_DIR(fastaprot)
        EUKULELE(MV_DIR.out.contigs_dir, EUKULELE_DB.out.db)

    emit:
       taxonomy_estimation = EUKULELE.out.taxonomy_extimation
       taxonomy_counts     = EUKULELE.out.taxonomy_counts
       diamond             = EUKULELE.out.diamond
       
       versions            = EUKULELE.out.versions

}
