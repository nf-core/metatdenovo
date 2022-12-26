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

// Validate parameters for assembler:
RNASPADES = 'rnaspades'
MEGAHIT   = 'megahit'

// validate parameters for eukulele database:
EUKULELE_DB_PHYLODB     = 'phylodb'
EUKULELE_DB_MMETSP      = 'mmetsp'
EUKULELE_DB_EUKPROT     = 'eukprot'
EUKULELE_DB_EUKZOO      = 'eukzoo'

def valid_params = [
    orf_caller      : [ORF_CALLER_PRODIGAL, ORF_CALLER_PROKKA, ORF_CALLER_TRANSDECODER],
    assembler       : [RNASPADES, MEGAHIT],
    eukulele_db     : [EUKULELE_DB_PHYLODB, EUKULELE_DB_MMETSP, EUKULELE_DB_EUKPROT, EUKULELE_DB_EUKZOO]
]

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// If the user supplied hmm files, we will run hmmsearch and then rank the results.
// Create a channel for hmm files.
if ( params.hmmsearch ) {
    if ( params.hmmdir ) {
        Channel
            .fromPath(params.hmmdir + params.hmmpattern)
            .set { ch_hmmrs }
    } else if ( params.hmmfiles ) {
        Channel
            .fromPath(params.hmmfiles)
            .set { ch_hmmrs }
    } else {
        // Warn of missing params
    }
}

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
include { HMMRANK                          } from '../modules/local/hmmrank.nf'
include { UNPIGZ as UNPIGZ_FASTA_PROTEIN   } from '../modules/local/unpigz.nf'
include { UNPIGZ as UNPIGZ_EUKULELE        } from '../modules/local/unpigz.nf'
include { UNPIGZ as UNPIGZ_EGGNOG          } from '../modules/local/unpigz.nf'
include { UNPIGZ as UNPIGZ_CONTIGS         } from '../modules/local/unpigz.nf'
include { COLLECT_FEATURECOUNTS            } from '../modules/local/collect_featurecounts.nf'
include { COLLECT_FEATURECOUNTS_EUK        } from '../modules/local/collect_featurecounts_euk.nf'
include { COLLECT_STATS                    } from '../modules/local/collect_stats.nf'
include { COLLECT_STATS_NOTRIM             } from '../modules/local/collect_stats_notrim.nf'

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

//
// SUBWORKFLOW: Consisting of local modules
//

include { EGGNOG        } from '../subworkflows/local/eggnog'
include { SUB_EUKULELE  } from '../subworkflows/local/eukulele'
include { EUKULELE_SRUN } from '../subworkflows/local/eukulele_second_run'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { BBMAP_BBDUK                                } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_INDEX                                } from '../modules/nf-core/bbmap/index/main'
include { BBMAP_ALIGN                                } from '../modules/nf-core/bbmap/align/main'
include { SEQTK_MERGEPE                              } from '../modules/nf-core/seqtk/mergepe/main'
include { BAM_SORT_SAMTOOLS                          } from '../subworkflows/nf-core/bam_sort_samtools/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/subread/featurecounts/main'
include { PRODIGAL                                   } from '../modules/nf-core/prodigal/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { HMMER_HMMSEARCH as HMMSEARCH               } from '../modules/nf-core/hmmer/hmmsearch/main.nf'
include { SPADES                                     } from '../modules/nf-core/spades/main'

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
        ch_versions = ch_versions.mix(DIGINORM.out.versions)
        ch_pe_reads_to_assembly = DIGINORM.out.pairs
        ch_se_reads_to_assembly = DIGINORM.out.singles
    } else {
        ch_pe_reads_to_assembly = SEQTK_MERGEPE.out.reads.map { meta, fastq -> fastq }
        ch_se_reads_to_assembly = []
    }

    ch_pacbio = []
    ch_nanopore = []
    ch_hmm = [] 
    ch_spades = SEQTK_MERGEPE.out.reads.map { [ [ id: 'all_samples' ], it[1],  [], [] ] } 
    

    //
    // MODULE: Run Megahit or RNAspades on all interleaved fastq files
    //
    if ( params.assembler == RNASPADES ) {
        ch_spades = FASTQC_TRIMGALORE.out.reads.map { meta, fastq -> [ [ id: 'all_samples' ], fastq, [], [] ] }
        SPADES( ch_spades, [] )
        ch_assembly_contigs = SPADES.out.transcripts.map { it[1] }
        ch_assembly_contigs.view()
        ch_versions = ch_versions.mix(SPADES.out.versions)
    } 
    if ( params.assembler == MEGAHIT ) {
    MEGAHIT_INTERLEAVED(
        ch_pe_reads_to_assembly.collect(),
        ch_se_reads_to_assembly.collect(),
        'all_samples'
    )
    ch_assembly_contigs = MEGAHIT_INTERLEAVED.out.contigs
    ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)
    }

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
        PROKKA_CAT(ch_assembly_contigs)
        ch_versions = ch_versions.mix(PROKKA_CAT.out.versions)
        ch_gff      = PROKKA_CAT.out.gff.map { it[1] }
        ch_hmm_aa   = PROKKA_CAT.out.faa
        ch_aa       = PROKKA_CAT.out.faa.map { it[1] }
    }

    //
    // MODULE: Call Prodigal
    //

    ch_prodigal = Channel.empty()
    if( params.orf_caller == ORF_CALLER_PRODIGAL ) {

        UNPIGZ_CONTIGS(ch_assembly_contigs)
        ch_versions = ch_versions.mix(UNPIGZ_CONTIGS.out.versions)
        PRODIGAL(
            UNPIGZ_CONTIGS.out.unzipped.collect { [ [ id: 'all_samples' ], it ] },
            'gff'
        )
        ch_gff          = PRODIGAL.out.gene_annotations.map { it[1] }
        ch_hmm_aa       = PRODIGAL.out.amino_acid_fasta
        ch_aa           = PRODIGAL.out.amino_acid_fasta
        ch_prodigal_fna = PRODIGAL.out.nucleotide_fasta
        ch_eukulele     = PRODIGAL.out.amino_acid_fasta
        ch_versions     = ch_versions.mix(PRODIGAL.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // SUBWORKFLOW: run TRANSDECODER on UNPIGZ output. Orf caller alternative for eukaryotes.
    //

    ch_transdecoder_longorf = Channel.empty()
    if( params.orf_caller == ORF_CALLER_TRANSDECODER ) {
        UNPIGZ_CONTIGS(ch_assembly_contigs)
        TRANSDECODER(
            UNPIGZ_CONTIGS.out.unzipped.collect { [ [ id: 'all_samples' ], it ] }
        )
        ch_gff      = TRANSDECODER.out.gff.map { it[1] }
        ch_hmm_aa   = TRANSDECODER.out.pep
        ch_aa       = TRANSDECODER.out.pep
        ch_gff      = TRANSDECODER.out.gff.map { it[1] }
        ch_eukulele = TRANSDECODER.out.pep
        ch_versions = ch_versions.mix(TRANSDECODER.out.versions)
    }

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //

    if (params.eggnog) {
        if ( params.orf_caller == ORF_CALLER_PROKKA ) {
            UNPIGZ_EGGNOG(ch_aa)
            EGGNOG(UNPIGZ_EGGNOG.out.unzipped.collect { [ [ id: 'all_samples' ], it ] } )
            ch_versions = ch_versions.mix(EGGNOG.out.versions)
        } else
        EGGNOG(ch_aa)
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
    }

    //
    // MODULE: Hmmsearch on orf caller output
    //
    if( params.hmmsearch) {
        ch_hmmstage = ch_hmmrs.combine(ch_hmm_aa.map { it[1] } )
            .map { [ [id: it[0].baseName ], it[0], it[1], true, true, false ] }
            .set { ch_hmmdir }
        HMMSEARCH( ch_hmmdir )
        HMMRANK( HMMSEARCH.out.target_summary.collect() { it[1] } )
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
    // MODULE: Collect featurecounts output counts in one table
    //

    if ( params.orf_caller == ORF_CALLER_PROKKA) {
        COLLECT_FEATURECOUNTS ( FEATURECOUNTS_CDS.out.counts.collect() { it[1] })
        ch_cds_counts = COLLECT_FEATURECOUNTS.out.counts
        ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    } else if ( params.orf_caller == ORF_CALLER_PRODIGAL) {
        COLLECT_FEATURECOUNTS ( FEATURECOUNTS_CDS.out.counts.collect() { it[1] })
        ch_cds_counts = COLLECT_FEATURECOUNTS.out.counts
        ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    } else if ( params.orf_caller == ORF_CALLER_TRANSDECODER) {
        COLLECT_FEATURECOUNTS_EUK ( FEATURECOUNTS_CDS.out.counts.collect() { it[1] })
        ch_cds_counts = COLLECT_FEATURECOUNTS_EUK.out.counts
        ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTS_EUK.out.versions)
    }
    ch_fcs = Channel.empty()
    ch_fcs = ch_fcs.mix(
        ch_cds_counts).collect()
    
    //
    // MODULE: Collect statistics from mapping analysis
    //

    if ( ! params.skip_trimming) {
        COLLECT_STATS (
            FASTQC_TRIMGALORE.out.trim_log.map { meta, fastq -> meta.id }.collect(),
            FASTQC_TRIMGALORE.out.trim_log.map { meta, fastq -> fastq[0] }.collect(),
            BAM_SORT_SAMTOOLS.out.idxstats.collect()  { it[1] },
            ch_fcs,
            ch_bbduk_logs.collect()
        )
        ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)
    } else {
        COLLECT_STATS_NOTRIM (
            FASTQC_TRIMGALORE.out.fastqc_html.map { meta, fastq -> meta.id }.collect(),
            BAM_SORT_SAMTOOLS.out.idxstats.collect()  { it[1] },
            ch_fcs,
            ch_bbduk_logs.collect()
        )
        ch_versions     = ch_versions.mix(COLLECT_STATS_NOTRIM.out.versions)
    }

    //
    // SUBWORKFLOW: Eukulele
    //


    if( params.eukulele ){
        if ( params.orf_caller == ORF_CALLER_PROKKA ) {
            UNPIGZ_EUKULELE(ch_aa)
            SUB_EUKULELE(UNPIGZ_EUKULELE.out.unzipped.collect { [ [ id: 'all_samples' ], it ] } )
        } else
        SUB_EUKULELE(ch_eukulele)
        ch_versions = ch_versions.mix(SUB_EUKULELE.out.versions)
    }

    if( params.eukulele_multirun){
        EUKULELE_SRUN(ch_eukulele)
        ch_versions = ch_versions.mix(EUKULELE_SRUN.out.versions)
    }

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files,
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
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
*/
