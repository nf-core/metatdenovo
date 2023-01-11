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
EUKULELE_DB_GTDB        = 'gtdb'

def valid_params = [
    orf_caller      : [ ORF_CALLER_PRODIGAL, ORF_CALLER_PROKKA, ORF_CALLER_TRANSDECODER ],
    assembler       : [ RNASPADES, MEGAHIT ],
    eukulele_db     : [ EUKULELE_DB_PHYLODB, EUKULELE_DB_MMETSP, EUKULELE_DB_EUKPROT, EUKULELE_DB_EUKZOO, EUKULELE_DB_GTDB ]
]

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// If the user supplied hmm files, we will run hmmsearch and then rank the results.
// Create a channel for hmm files.
ch_hmmrs = Channel.empty()
if ( params.hmmdir ) {
    Channel
        .fromPath(params.hmmdir + params.hmmpattern)
        .set { ch_hmmrs }
} else if ( params.hmmfiles ) {
    Channel
        .of( params.hmmfiles.split(',') )
        .map { [ file(it) ] }
        .set { ch_hmmrs }
}

ch_eukulele_db = Channel.empty()
if ( !params.skip_eukulele ) {
    if ( params.eukulele_db) {
        Channel
            .of ( params.eukulele_db.split(',') )
            .set { ch_eukulele_db }
    } else {  exit 1, 'eukuelel database not specified!' 
    }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local
//
include { MEGAHIT_INTERLEAVED              } from '../modules/local/megahit/interleaved.nf'
include { UNPIGZ as UNPIGZ_FASTA_PROTEIN   } from '../modules/local/unpigz.nf'
include { UNPIGZ as UNPIGZ_EUKULELE        } from '../modules/local/unpigz.nf'
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
// SUBWORKFLOW
//

include { PROKKA_SUBSETS } from '../subworkflows/local/prokka_subsets'
include { TRANSDECODER   } from '../subworkflows/local/transdecoder'
include { EGGNOG        } from '../subworkflows/local/eggnog'
include { SUB_EUKULELE  } from '../subworkflows/local/eukulele'
include { HMMCLASSIFY   } from '../subworkflows/local/hmmclassify'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BBMAP_BBDUK                                } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_INDEX                                } from '../modules/nf-core/bbmap/index/main'
include { BBMAP_ALIGN                                } from '../modules/nf-core/bbmap/align/main'
include { SEQTK_MERGEPE                              } from '../modules/nf-core/seqtk/mergepe/main'
include { BAM_SORT_SAMTOOLS                          } from '../subworkflows/nf-core/bam_sort_samtools/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/subread/featurecounts/main'
include { PRODIGAL                                   } from '../modules/nf-core/prodigal/main'
include { SPADES                                     } from '../modules/nf-core/spades/main'
include { CAT_FASTQ 		                         } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'

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
    .reads
    .map {
        meta, fastq ->
            new_id = meta.id - ~/_T\d+/
            [ meta + [id: new_id], fastq ]
    }
    .groupTuple()
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // Branch FastQ channels if 'auto' specified to infer strandedness
    ch_cat_fastq
        .branch {
            meta, fastq ->
                auto_strand : meta.strandedness == 'auto'
                    return [ meta, fastq ]
                known_strand: meta.strandedness != 'auto'
                    return [ meta, fastq ]
        }
        .set { ch_strand_fastq }
    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
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
    // MODULE: Interleave sequences for assembly
    //
    ch_interleaved = Channel.empty()
    if ( ! params.assembly ) {
        SEQTK_MERGEPE(ch_clean_reads)
        ch_interleaved = SEQTK_MERGEPE.out.reads
        ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)
    }

    //
    // SUBWORKFLOW: Perform digital normalization
    //
    ch_pe_reads_to_assembly = Channel.empty()
    ch_se_reads_to_assembly = Channel.empty()
    if ( ! params.assembly ) {
        if ( params.diginorm ) {
            DIGINORM(ch_interleaved.collect { meta, fastq -> fastq }, [], 'all_samples')
            ch_versions = ch_versions.mix(DIGINORM.out.versions)
            ch_pe_reads_to_assembly = DIGINORM.out.pairs
            ch_se_reads_to_assembly = DIGINORM.out.singles
        } else {
            ch_pe_reads_to_assembly = ch_interleaved.map { meta, fastq -> fastq }
            ch_se_reads_to_assembly = []
        }
    }

    //
    // MODULE: Run Megahit or RNAspades on all interleaved fastq files
    //
    if ( params.assembly ) {
        Channel
            .value ( [ [ id: 'user_assembly' ], file(params.assembly) ] )
            .set { ch_assembly_contigs }
    } else if ( params.assembler == RNASPADES ) {
        // This doesn't work as we want, as it gets called once for each pair, see issue: https://github.com/LNUc-EEMiS/metatdenovo/issues/78
        ch_spades = FASTQC_TRIMGALORE.out.reads.map { meta, fastq -> [ [ id: 'spades' ], fastq, [], [] ] }
        SPADES( ch_spades, [], [] )
        ch_assembly_contigs = SPADES.out.transcripts.map { it[1] }
        ch_versions = ch_versions.mix(SPADES.out.versions)
    } else if ( params.assembler == MEGAHIT ) {
        MEGAHIT_INTERLEAVED(
            ch_pe_reads_to_assembly.collect(),
            ch_se_reads_to_assembly.collect(),
            'megahit_assembly'
        )
        MEGAHIT_INTERLEAVED.out.contigs
            .map { [ [ id: 'megahit' ], it ] }
            .set { ch_assembly_contigs }
        ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)
    }

    //
    // MODULE: Create a BBMap index
    //
    BBMAP_INDEX(ch_assembly_contigs.map { it[1] })
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
    // Call ORFs
    //

    ch_gff = Channel.empty()
    ch_aa  = Channel.empty()

    //
    // SUBWORKFLOW: Run PROKKA_SUBSETS on Megahit output, but split the fasta file in chunks of 10 MB, then concatenate and compress output.
    //

    if ( params.orf_caller == ORF_CALLER_PROKKA ) {
        PROKKA_SUBSETS(ch_assembly_contigs)
        ch_versions = ch_versions.mix(PROKKA_SUBSETS.out.versions)
        // DL: Isn't it clearer to leave the mapping to when the channel is used?
        ch_gff      = PROKKA_SUBSETS.out.gff.map { it[1] }
        ch_aa       = PROKKA_SUBSETS.out.faa
    }

    //
    // MODULE: Call Prodigal
    //

    if ( params.orf_caller == ORF_CALLER_PRODIGAL ) {

        UNPIGZ_CONTIGS(ch_assembly_contigs.map { it[1] })
        ch_versions = ch_versions.mix(UNPIGZ_CONTIGS.out.versions)

        PRODIGAL ( ch_assembly_contigs, 'gff' )
        ch_gff          = PRODIGAL.out.gene_annotations.map { it[1] }
        ch_aa           = PRODIGAL.out.amino_acid_fasta
        ch_versions     = ch_versions.mix(PRODIGAL.out.versions)
    }

    //
    // SUBWORKFLOW: run TRANSDECODER on UNPIGZ output. Orf caller alternative for eukaryotes.
    //

    if ( params.orf_caller == ORF_CALLER_TRANSDECODER ) {

        UNPIGZ_CONTIGS(ch_assembly_contigs.map { it[1] })
        ch_versions = ch_versions.mix(UNPIGZ_CONTIGS.out.versions)

        TRANSDECODER ( ch_assembly_contigs )

        ch_gff      = TRANSDECODER.out.gff.map { it[1] }
        ch_aa       = TRANSDECODER.out.pep
        ch_versions = ch_versions.mix(TRANSDECODER.out.versions)
    }

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //

    if (params.eggnog) {
        EGGNOG(ch_aa)
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
    }

    //
    // SUBWORKFLOW: classify ORFs with a set of hmm files
    //

    ch_hmmrs
        .combine(ch_aa)
        .map { [ [id: it[0].baseName ], it[0], it[2] ] }
        .set { ch_hmmclassify }
    HMMCLASSIFY ( ch_hmmclassify )
    ch_versions = ch_versions.mix(HMMCLASSIFY.out.versions)

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
    
    FEATURECOUNTS_CDS.out.counts
        .collect() { it[1] }
        .map { [ 'all_samples', it ] }
        .set { ch_collect_feature }

    COLLECT_FEATURECOUNTS ( ch_collect_feature )

    ch_fcs = Channel.empty()
    ch_fcs = COLLECT_FEATURECOUNTS.out.counts.collect()
    ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)

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

    if( !params.skip_eukulele){
        ch_eukulele_db
            .combine(ch_aa)
            .map{ [ [id:it[0] ], it[2] ] }
            .set { ch_eukulele }
        SUB_EUKULELE( ch_eukulele, ch_eukulele_db )
        ch_versions = ch_versions.mix(SUB_EUKULELE.out.versions)
    }

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //

    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetatdenovo.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
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
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
