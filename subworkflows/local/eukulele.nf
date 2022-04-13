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
        ch_versions = Channel.empty() 
        
        MV_DIR(fastaprot) 
        
        String directoryName = params.eukulele_dbpath
        File directory = new File(directoryName)
        
        if(! directory.exists()){
            directory.mkdir()
            EUKULELE_DB( )
            EUKULELE_DB.out.database.mklink('./eukulele')
            EUKULELE(MV_DIR.out.contigs_dir, EUKULELE_DB.out.database)
        } else {
            ch_database = Channel.fromPath(params.eukulele_dbpath)
            EUKULELE(MV_DIR.out.contigs_dir, ch_database)
        }
        

    emit:
       taxonomy_estimation = EUKULELE.out.taxonomy_extimation
       taxonomy_counts     = EUKULELE.out.taxonomy_counts
       diamond             = EUKULELE.out.diamond
       
       versions            = EUKULELE.out.versions

}
