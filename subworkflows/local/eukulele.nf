//
// Run EUKULELE on contigs that are moved into a directory
//

include { EUKULELE} from '../../modules/local/eukulele/main'
include { MV_DIR  } from '../../modules/local/mvdir.nf'
include { WGET_DB } from '../../modules/local/eukulele/download'

workflow SUB_EUKULELE {

    if( params.eukulele_dbpath == '~/eukulele') { 
        if( params.eukulele_db == 'phylodb') {
            input1 = "https://www.dropbox.com/s/dl/ahpqdkmjg6xekwg/phylodb_1.076.pep.fa.gz"
            input2 = "https://www.dropbox.com/s/dl/ak76j4po9tdcbof/phylodb_1.076.taxonomy.txt"
        } else {
            input1 = "https://www.dropbox.com/s/4qilp6se3qf47uq/reference-pep.fa"
            input2 = "https://www.dropbox.com/s/w1uv1j9qr63z3ac/taxonomy_table.txt"
        }
    }

    take:
        contigs
        db

    main:
        ch_versions = Channel.empty()
        ch_eukulele_pathdb = Channel.empty()
        
        if( params.eukulele_dbpath == '~/eukulele'){ 
            WGET_DB         ( input1, input2 )
        }
        
        MV_DIR          ( contigs )
        
        if( params.eukulele_dbpath == '~/eukulele') {
            EUKULELE        ( MV_DIR.out.contigs_dir, WGET_DB.out.database )
        } 
        else {
            EUKULELE (MV_DIR.out.contigs_dir, db)
        }
        ch_versions = ch_versions.mix(EUKULELE.out.versions)
    
    emit:
    taxonomy_extimation = EUKULELE.out.taxonomy_extimation
    taxonomy_counts     = EUKULELE.out.taxonomy_counts
    diamond             = EUKULELE.out.diamond
        
       // versions            = ch_versions
}
