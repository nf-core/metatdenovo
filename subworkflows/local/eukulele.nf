//
// Run EUKULELE on protein fasta from orf_caller output
//

include { EUKULELE      } from '../../modules/local/eukulele/main'
include { MV_DIR        } from '../../modules/local/mvdir'
include { WGET_DB       } from '../../modules/local/eukulele/download'
include { TABLE_PHYLODB } from '../../modules/local/eukulele/create_table_phylodb'
include { TABLE_MMETSP  } from '../../modules/local/eukulele/create_table_mmetsp'

workflow SUB_EUKULELE {

    if( params.eukulele_db == 'phylodb') {
        input1= file("https://www.dropbox.com/s/ahpqdkmjg6xekwg/phylodb_1.076.pep.fa.gz?dl=1")
        input2= file("https://www.dropbox.com/s/ak76j4po9tdcbof/phylodb_1.076.taxonomy.txt.gz?dl=1")
    } else {
        input1= file("https://www.dropbox.com/s/4qilp6se3qf47uq/reference-pep.fa?dl=1")
        input2= file("https://www.dropbox.com/s/w1uv1j9qr63z3ac/taxonomy_table.txt?dl=1")
    }

    take:
        fastaprot

    main:
        
        WGET_DB(input1,input2)
        if( params.eukulele_db == 'phylodb' ){
            TABLE_PHYLODB(WGET_DB.out.database, WGET_DB.out.taxonomy)
            MV_DIR(fastaprot)
            EUKULELE(MV_DIR.out.contigs_dir, TABLE_PHYLODB.out.db)
        } else {
            TABLE_MMETSP(WGET_DB.out.database, WGET_DB.out.taxonomy)
            MV_DIR(fastaprot)
            EUKULELE(MV_DIR.out.contigs_dir, TABLE_MMETSP.out.db)
        }

    emit:
       taxonomy_estimation = EUKULELE.out.taxonomy_extimation
       taxonomy_counts     = EUKULELE.out.taxonomy_counts
       diamond             = EUKULELE.out.diamond
       
       versions            = EUKULELE.out.versions

}
