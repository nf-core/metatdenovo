/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMetatdenovo.initialise(params, log)

// set an empty multiqc channel
ch_multiqc_files = Channel.empty()

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

// Create a channel for EUKulele either with a named database or not. The latter means a user-provided database in a directory.
ch_eukulele_db = Channel.empty()
if ( !params.skip_eukulele ) {
    if ( params.eukulele_db ) {
        Channel
            .of ( params.eukulele_db.split(',') )
            .map { [ it, file(params.eukulele_dbpath) ] }
            .set { ch_eukulele_db }
    } else {
        Channel.fromPath(params.eukulele_dbpath)
            .map { [ [], it ] }
            .set { ch_eukulele_db }
    }
}

// Create a channel for CAT that contains the path for the database
if(params.cat_db){
    ch_cat_db_file = Channel
        .value(file( "${params.cat_db}" ))
} else {
    ch_cat_db_file = Channel.empty()
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
include { WRITESPADESYAML                  } from '../modules/local/writespadesyaml'
include { MEGAHIT_INTERLEAVED              } from '../modules/local/megahit/interleaved'
include { COLLECT_FEATURECOUNTS            } from '../modules/local/collect_featurecounts'
include { COLLECT_STATS                    } from '../modules/local/collect_stats'
include { CAT_DB                           } from '../modules/local/cat/cat_db'
include { CAT_DB_GENERATE                  } from '../modules/local/cat/cat_db_generate'
include { CAT_CONTIGS                      } from '../modules/local/cat/cat_contigs'
include { CAT_SUMMARY                      } from "../modules/local/cat/cat_summary"
include { FORMATSPADES                     } from '../modules/local/formatspades'
include { UNPIGZ as UNPIGZ_CONTIGS         } from '../modules/local/unpigz'
include { UNPIGZ as UNPIGZ_GFF             } from '../modules/local/unpigz'
include { MERGE_TABLES                     } from '../modules/local/merge_summary_tables'
include { TRANSRATE                        } from '../modules/local/transrate'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// SUBWORKFLOW: Consisting of local modules
//
include { EGGNOG            } from '../subworkflows/local/eggnog'
include { SUB_EUKULELE      } from '../subworkflows/local/eukulele'
include { HMMCLASSIFY       } from '../subworkflows/local/hmmclassify'
include { PROKKA_SUBSETS    } from '../subworkflows/local/prokka_subsets'
include { TRANSDECODER      } from '../subworkflows/local/transdecoder'
include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore'
include { PRODIGAL          } from '../subworkflows/local/prodigal'
include { KOFAMSCAN         } from '../subworkflows/local/kofamscan'

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
include { BBMAP_BBNORM                               } from '../modules/nf-core/bbmap/bbnorm/main'
include { SEQTK_MERGEPE                              } from '../modules/nf-core/seqtk/mergepe/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/subread/featurecounts/main'
include { SPADES                                     } from '../modules/nf-core/spades/main'
include { SEQTK_SEQ as SEQTK_SEQ_CONTIG_FILTER       } from '../modules/nf-core/seqtk/seq/main'
include { CAT_FASTQ            	                     } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { BAM_SORT_STATS_SAMTOOLS                    } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

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
        file(params.input)
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
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

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

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    FASTQC_TRIMGALORE.out.trim_log.collect { it[1] }
    ch_collect_stats = ch_cat_fastq.collect { it[0].id }.map { [ [ id:"${params.assembler}.${params.orf_caller}" ], it ] }
    if ( params.skip_trimming ) {
        ch_collect_stats
            .map { [ it[0], it[1], [] ] }
            .set { ch_collect_stats }
    } else {
        if ( params.se_reads ) {
            ch_collect_stats
                .combine(FASTQC_TRIMGALORE.out.trim_log.collect { it[1] }.map { [ it ] })
                .set { ch_collect_stats }
        } else {
            ch_collect_stats
                .combine(FASTQC_TRIMGALORE.out.trim_log.collect { it[1][0] }.map { [ it ] })
                .set { ch_collect_stats }
        }
    }

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, params.sequence_filter )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { it[1] }.map { [ it ] }
        ch_versions   = ch_versions.mix(BBMAP_BBDUK.out.versions)
        ch_collect_stats
            .combine(ch_bbduk_logs)
            .set {ch_collect_stats}
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{it[1]}.ifEmpty([]))
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = Channel.empty()
        ch_collect_stats
            .map { [ it[0], it[1], it[2], [] ] }
            .set { ch_collect_stats }
    }

    //
    // MODULE: Interleave sequences for assembly
    //
    // DL & DDL: We can probably not deal with single end input
    ch_interleaved = Channel.empty()
    if ( ! params.assembly ) {
        if ( params.se_reads) {
            ch_single_end = ch_clean_reads
        } else {
            SEQTK_MERGEPE(ch_clean_reads)
            ch_interleaved = SEQTK_MERGEPE.out.reads
            ch_versions    = ch_versions.mix(SEQTK_MERGEPE.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Perform digital normalization. There are two options: khmer or BBnorm. The latter is faster.
    //
    if ( ! params.assembly ) {
        if ( params.se_reads ) {
            if ( params.bbnorm ) {
                BBMAP_BBNORM(ch_single_end.collect { meta, fastq -> fastq }.map {[ [id:'all_samples', single_end:true], it ] } )
                ch_se_reads_to_assembly = BBMAP_BBNORM.out.fastq.map { it[1] }
                ch_pe_reads_to_assembly = Channel.empty()
            } else {
                ch_se_reads_to_assembly = ch_single_end.map { meta, fastq -> fastq }
                ch_pe_reads_to_assembly = Channel.empty()
            }
        }
        else if ( params.bbnorm ) {
            BBMAP_BBNORM(ch_interleaved.collect { meta, fastq -> fastq }.map {[ [id:'all_samples', single_end:true], it ] } )
            ch_pe_reads_to_assembly = BBMAP_BBNORM.out.fastq.map { it[1] }
            ch_se_reads_to_assembly = Channel.empty()
        } else {
            ch_pe_reads_to_assembly = ch_interleaved.map { meta, fastq -> fastq }
            ch_se_reads_to_assembly = Channel.empty()
        }
    }

    //
    // MODULE: Run Megahit or RNAspades on all interleaved fastq files
    //
    if ( params.assembly ) {
        Channel
            .value ( [ [ id: 'user_assembly' ], file(params.assembly) ] )
            .set { ch_assembly_contigs }
    } else if ( params.assembler == 'rnaspades' ) {
        // 1. Write a yaml file for Spades
        WRITESPADESYAML (
            ch_pe_reads_to_assembly.collect().ifEmpty([]),
            ch_se_reads_to_assembly.collect().ifEmpty([])
        )
        // 2. Call the module with a channel with all fastq files plus the yaml
        Channel.empty()
            .mix(ch_pe_reads_to_assembly)
            .mix(ch_se_reads_to_assembly)
            .collect()
            .map { [ [ id:'rnaspades' ], it, [], [] ] }
            .set { ch_spades }
        SPADES (
            ch_spades,
            WRITESPADESYAML.out.yaml,
            []
        )
        ch_assembly = SPADES.out.transcripts
        ch_versions = ch_versions.mix(SPADES.out.versions)
        FORMATSPADES( ch_assembly )
        ch_assembly_contigs = FORMATSPADES.out.assembly
    } else if ( params.assembler == 'megahit' ) {
        MEGAHIT_INTERLEAVED(
            ch_pe_reads_to_assembly.collect().ifEmpty([]),
            ch_se_reads_to_assembly.collect().ifEmpty([]),
            'megahit_assembly'
        )
        MEGAHIT_INTERLEAVED.out.contigs
            .map { [ [ id: 'megahit' ], it ] }
            .set { ch_assembly_contigs }
        ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)
    } else { exit 1, 'Assembler not specified!' }

    // If the user asked for length filtering, perform that with SEQTK_SEQ (the actual length parameter is used in modules.config)
    if ( params.min_contig_length > 0 ) {
        SEQTK_SEQ_CONTIG_FILTER ( ch_assembly_contigs )
        ch_assembly_contigs = SEQTK_SEQ_CONTIG_FILTER.out.fastx
    }

    //
    // Call ORFs
    //
    ch_gff = Channel.empty()
    ch_aa  = Channel.empty()

    //
    // SUBWORKFLOW: Run PROKKA_SUBSETS on assmebly output, but split the fasta file in chunks of 10 MB, then concatenate and compress output.
    //
    if ( params.orf_caller == 'prokka' ) {
        PROKKA_SUBSETS(ch_assembly_contigs, params.prokka_batchsize)
        UNPIGZ_GFF(PROKKA_SUBSETS.out.gff.map { [ [id: "${params.orf_caller}.${it[0].id}"], it[1] ] })
        ch_versions      = ch_versions.mix(PROKKA_SUBSETS.out.versions)
        ch_gff           = UNPIGZ_GFF.out.unzipped
        ch_aa            = PROKKA_SUBSETS.out.faa
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA_SUBSETS.out.prokka_log.collect{it[1]}.ifEmpty([]))
    }

    //
    // MODULE: Run PRODIGAL on assembly output.
    //
    if ( params.orf_caller == 'prodigal' ) {
        PRODIGAL( ch_assembly_contigs.map { [ [id: "${params.assembler}.${params.orf_caller}"], it[1] ] } )
        UNPIGZ_GFF(PRODIGAL.out.gff.map { [ [id: "${it[0].id}.${params.orf_caller}"], it[1] ] })
        ch_gff          = UNPIGZ_GFF.out.unzipped
        ch_aa           = PRODIGAL.out.faa
        ch_versions     = ch_versions.mix(PRODIGAL.out.versions)
    }

    //
    // SUBWORKFLOW: run TRANSDECODER. Orf caller alternative for eukaryotes.
    //
    if ( params.orf_caller == 'transdecoder' ) {
        TRANSDECODER ( ch_assembly_contigs.map { [ [id: "transdecoder.${it[0].id}" ], it[1] ] } )
        ch_gff      = TRANSDECODER.out.gff
        ch_aa       = TRANSDECODER.out.pep
        ch_versions = ch_versions.mix(TRANSDECODER.out.versions)
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
    // SUBWORKFLOW: classify ORFs with a set of hmm files
    //
    ch_hmmrs
        .combine(ch_aa)
        .map { [ [ id: "${params.assembler}.${params.orf_caller}" ], it[0], it[2] ] }
        .set { ch_hmmclassify }
    HMMCLASSIFY ( ch_hmmclassify )
    ch_versions = ch_versions.mix(HMMCLASSIFY.out.versions)

    //
    // MODULE: FeatureCounts. Create a table for each samples that provides raw counts as result of the alignment.
    //

    BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, ch_assembly_contigs )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // if ( orf_caller ==
    BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(ch_gff.map { it[1] } )
        .set { ch_featurecounts }

    ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it[1]}.map { [ it ] })
        .set { ch_collect_stats }

    FEATURECOUNTS_CDS ( ch_featurecounts)
    ch_versions       = ch_versions.mix(FEATURECOUNTS_CDS.out.versions)

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    FEATURECOUNTS_CDS.out.counts
        .collect() { it[1] }
        .map { [ [ id:"${params.assembler}.${params.orf_caller}" ], it ] }
        .set { ch_collect_feature }

    COLLECT_FEATURECOUNTS ( ch_collect_feature )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { it[1]}.map { [ it ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { it[1]}
    ch_collect_stats
        .combine(ch_fcs_for_stats)
        .set { ch_collect_stats }

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //
    if ( ! params.skip_eggnog ) {
        File directory       = new File(params.eggnog_dbpath)
        if ( ! directory.exists() ) { directory.mkdir() }
        EGGNOG(params.eggnog_dbpath, ch_aa, ch_fcs_for_summary )
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
        ch_merge_tables = EGGNOG.out.sumtable
    } else {
        ch_aa
            .map { [ it[0], [] ] }
            .set { ch_merge_tables }
    }


    //
    // SUBWORKFLOW: run kofamscan on the ORF-called amino acid sequences
    //
    if( !params.skip_kofamscan ) {
        ch_aa
            .map { [ it[0], it[1] ] }
            .set { ch_kofamscan }
        KOFAMSCAN( ch_kofamscan, params.kofam_dir, ch_fcs_for_summary)
        ch_versions = ch_versions.mix(KOFAMSCAN.out.versions)
        ch_kofamscan_summary = KOFAMSCAN.out.kofamscan_summary.collect().map { it[1] }
        ch_merge_tables
            .combine( ch_kofamscan_summary )
            .set { ch_merge_tables }
    } else {
        ch_merge_tables
            .map { [ it[0], it[1], [] ] }
            .set { ch_merge_tables }
    }

    // set up contig channel to use in CAT and TransRate
    UNPIGZ_CONTIGS(ch_assembly_contigs)
    ch_unzipped_contigs = UNPIGZ_CONTIGS.out.unzipped
    ch_versions = ch_versions.mix(UNPIGZ_CONTIGS.out.versions)

    //
    // CAT: Bin Annotation Tool (BAT) are pipelines for the taxonomic classification of long DNA sequences and metagenome assembled genomes (MAGs/bins)
    //
    ch_cat_db = Channel.empty()
    if (params.cat) {
        if (params.cat_db){
            CAT_DB ( ch_cat_db_file )
            ch_cat_db = CAT_DB.out.db
        } else if (params.cat_db_generate){
            CAT_DB_GENERATE ()
            ch_cat_db = CAT_DB_GENERATE.out.db
        }
        CAT_CONTIGS (
            ch_unzipped_contigs,
            ch_cat_db
        )
        CAT_SUMMARY(
            CAT_CONTIGS.out.tax_classification.collect()
        )
        ch_versions = ch_versions.mix(CAT_CONTIGS.out.versions)
        ch_versions = ch_versions.mix(CAT_SUMMARY.out.versions)
    }

    //
    // MODULE: Use TransRate to judge assembly quality, piped into MultiQC
    //
    TRANSRATE(ch_unzipped_contigs)
    ch_versions = ch_versions.mix(TRANSRATE.out.versions)

    //
    // SUBWORKFLOW: Eukulele
    //
    if( !params.skip_eukulele){
        File directory = new File(params.eukulele_dbpath)
        if ( ! directory.exists() ) { directory.mkdir() }
        ch_directory = Channel.fromPath( directory )
            ch_aa
                .map {[ [ id:"${it[0].id}" ], it[1] ] }
                .combine( ch_eukulele_db )
                .set { ch_eukulele }
            SUB_EUKULELE( ch_eukulele, ch_fcs_for_summary )
            ch_taxonomy_summary = SUB_EUKULELE.out.taxonomy_summary.collect().map { it[1] }
            ch_versions = ch_versions.mix(SUB_EUKULELE.out.versions)
            ch_merge_tables
                .combine( ch_taxonomy_summary )
                .set { ch_merge_tables }
    } else {
        ch_merge_tables
            .map { [ it[0], it[1], it[2], [] ] }
            .set { ch_merge_tables }
    }

    //
    // MODULE: Collect statistics from mapping analysis
    //
    if( !params.skip_eggnog  || !params.skip_eukulele || !params.skip_kofamscan) {
        MERGE_TABLES ( ch_merge_tables )
        ch_collect_stats
            .combine(MERGE_TABLES.out.merged_table.collect{ it[1]}.map { [ it ] })
            .set { ch_collect_stats }
    } else {
        ch_collect_stats
            .map { [ it[0], it[1], it[2], it[3], it[4], it[5], [] ] }
            .set { ch_collect_stats }
    }

    COLLECT_STATS(ch_collect_stats)
    ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetatdenovo.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(TRANSRATE.out.assembly_qc.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS_CDS.out.summary.collect{it[1]}.ifEmpty([]))


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
    NfcoreTemplate.dump_parameters(workflow, params)
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
