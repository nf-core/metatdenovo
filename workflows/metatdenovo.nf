/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetatdenovo.initialise(params, log)

// Validate parameters for orf_caller:
ORF_CALLER_PRODIGAL     = 'prodigal'
ORF_CALLER_PROKKA       = 'prokka'
ORF_CALLER_TRANSDECODER = 'transdecoder'

// validate parameters for eukulele database:
EUKULELE_DB_PHYLODB     = 'phylodb'
EUKULELE_DB_MMETSP      = 'mmetsp'
EUKULELE_DB_EUKPROT     = 'eukprot'
EUKULELE_DB_EUKZOO      = 'eukzoo'

def valid_params = [
    orf_caller  : [ORF_CALLER_PRODIGAL, ORF_CALLER_PROKKA, ORF_CALLER_TRANSDECODER],
    eukulele_db : [EUKULELE_DB_PHYLODB, EUKULELE_DB_MMETSP, EUKULELE_DB_EUKPROT, EUKULELE_DB_EUKZOO ]
]

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULE: local
//
include { MEGAHIT_INTERLEAVED              } from '../modules/local/megahit/interleaved.nf'
include { UNPIGZ as UNPIGZ_MEGAHIT_CONTIGS } from '../modules/local/unpigz.nf'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// SUBWORKFLOW: Adapted from rnaseq!
//

include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore'

//
// SUBWORKFLOW: Perform digital normalization
//

include { DIGINORM } from '../subworkflows/local/diginorm'

//
// SUBWORKFLOW: Consisting of nf-core/modules
//

include { PROKKA_CAT   } from '../subworkflows/local/prokka_cat'
include { TRANSDECODER } from '../subworkflows/local/transdecoder'
include { SUB_EUKULELE } from '../subworkflows/local/eukulele'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                     } from '../modules/nf-core/modules/fastqc/main'
include { BBMAP_BBDUK                                } from '../modules/nf-core/modules/bbmap/bbduk/main'
include { BBMAP_INDEX                                } from '../modules/nf-core/modules/bbmap/index/main'
include { BBMAP_ALIGN                                } from '../modules/nf-core/modules/bbmap/align/main'
include { SEQTK_MERGEPE                              } from '../modules/nf-core/modules/seqtk/mergepe/main'
include { BAM_SORT_SAMTOOLS                          } from '../subworkflows/nf-core/bam_sort_samtools/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/modules/subread/featurecounts/main'
include { PRODIGAL                                   } from '../modules/nf-core/modules/prodigal/main'
include { MULTIQC                                    } from '../modules/nf-core/modules/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METATDENOVO {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, params.sequence_filter )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.map { it[1] }
        ch_versions   = ch_versions.mix(BBMAP_BBDUK.out.versions)
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = []
    }

    //
    // MODULE: Interleave sequences
    //
    SEQTK_MERGEPE(ch_clean_reads)
    ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)

    //
    // SUBWORKFLOW: Perform digital normalization
    //
    ch_reads_to_assembly = Channel.empty()
    if ( params.diginorm ) {
        DIGINORM(SEQTK_MERGEPE.out.reads.collect { meta, fastq -> fastq }, [], 'all_samples')
        ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)
        ch_pe_reads_to_assembly = DIGINORM.out.pairs
        ch_se_reads_to_assembly = DIGINORM.out.singles
    } else {
        ch_pe_reads_to_assembly = SEQTK_MERGEPE.out.reads.map { meta, fastq -> fastq }
        ch_se_reads_to_assembly = []
    }

    //
    // MODULE: Run Megahit on all interleaved fastq files
    //
    MEGAHIT_INTERLEAVED(
        ch_pe_reads_to_assembly.collect(),
        ch_se_reads_to_assembly.collect(),
        'all_samples'
    )
    ch_assembly_contigs = MEGAHIT_INTERLEAVED.out.contigs
    ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)

    //
    // MODULE: Create a BBMap index
    //
    BBMAP_INDEX(ch_assembly_contigs)
    ch_versions   = ch_versions.mix(BBMAP_INDEX.out.versions)

    //
    // MODULE: Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, BBMAP_INDEX.out.index )
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    //
    // SUBWORKFLOW: sort bam file
    //
    BAM_SORT_SAMTOOLS ( BBMAP_ALIGN.out.bam )
    ch_versions = ch_versions.mix(BAM_SORT_SAMTOOLS.out.versions)

    //
    // SUBWORKFLOW: Run PROKKA on Megahit output, but split the fasta file in chunks of 10 MB, then concatenate and compress output.
    //
    if (params.orf_caller == ORF_CALLER_PROKKA) {
        PROKKA_CAT(MEGAHIT_INTERLEAVED.out.contigs)
        ch_versions = ch_versions.mix(PROKKA_CAT.out.versions)
        ch_gff      = PROKKA_CAT.out.gff
        ch_eukulele = PROKKA_CAT.out.faa
    }

    //
    // MODULE: Call Prodigal
    //
    UNPIGZ_MEGAHIT_CONTIGS(ch_assembly_contigs)
    ch_versions = ch_versions.mix(UNPIGZ_MEGAHIT_CONTIGS.out.versions)

    ch_prodigal = Channel.empty()
    if( params.orf_caller == ORF_CALLER_PRODIGAL ) {
        PRODIGAL(
            UNPIGZ_MEGAHIT_CONTIGS.out.unzipped.collect { [ [ id: 'all_samples' ], it ] },
            'gff'
        )
        ch_gff          = PRODIGAL.out.gene_annotations.map { it[1] }
        ch_prodigal_aa  = PRODIGAL.out.amino_acid_fasta
        ch_prodigal_fna = PRODIGAL.out.nucleotide_fasta
        ch_eukulele     = PRODIGAL.out.amino_acid_fasta
        ch_versions     = ch_versions.mix(PRODIGAL.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // SUBWORKFLOW: run TRANSDECODER on UNPIGZ_MEGAHIT output. Orf caller alternative for eukaryotes.
    //

    ch_transdecoder_longorf = Channel.empty()
    if( params.orf_caller == ORF_CALLER_TRANSDECODER ) {
        TRANSDECODER(
            UNPIGZ_MEGAHIT_CONTIGS.out.unzipped.collect { [ [ id: 'all_samples' ], it ] }
        )
        ch_gff      = TRANSDECODER.out.gff.map { it[1] }
        ch_eukulele = TRANSDECODER.out.pep
    }

    //
    // MODULE: FeatureCounts
    //

    BAM_SORT_SAMTOOLS.out.bam
        .combine(ch_gff)
        .set { ch_featurecounts }

    FEATURECOUNTS_CDS ( ch_featurecounts)
    ch_versions       = ch_versions.mix(FEATURECOUNTS_CDS.out.versions)
    
    //
    // SUBWORKFLOW: Eukulele
    //
    
    if (params.eukulele_pathdb) {
        ch_eukulele_pathdb = Channel.fromPath(params.eukulele_pathdb)
    }
    else {
        ch_eukulele_pathdb = Channel.fromPath('.')
    }
    
    ch_unpigz_contigs = UNPIGZ_MEGAHIT_CONTIGS.out.unzipped.collect { [ [ id: 'all_samples' ], it ] }
    
    if( !params.skip_eukulele){
            SUB_EUKULELE(ch_eukulele, ch_eukulele_pathdb)
    }
    
    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    
    // Make sure we integrate FASTQC output from FASTQC_TRIMGALORE here!!!
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
