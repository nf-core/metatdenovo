/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local
//
include { COLLECT_FEATURECOUNTSSUMMARY       } from '../modules/local/collect/featurecountssummary/'
include { COLLECT_LOCUSCONSOLIDATE           } from '../modules/local/collect/locusconsolidate/'
include { COLLECT_PROTEINCONSOLIDATE         } from '../modules/local/collect/proteinconsolidate/'
include { FORMAT_CLUSTERREPS                 } from '../modules/local/format/clusterreps/'
include { FORMAT_GFF2BED                     } from '../modules/local/format/gff2bed/'
include { FORMAT_LOCUSCONSOLIDATE            } from '../modules/local/format/locusconsolidate/'
include { FORMAT_LOCUSFAA                    } from '../modules/local/format/locusfaa/'
include { FORMATSPADES                       } from '../modules/local/format/spades/'
include { MERGE_TABLES                       } from '../modules/local/merge/summary/'
include { FORMAT_DIAMOND_TAX_RANKLIST        } from '../modules/local/diamond/format_tax/ranklist/'
include { FORMAT_DIAMOND_TAX_TAXDUMP         } from '../modules/local/diamond/format_tax/taxdump/'
include { SUMTAXONOMY as SUM_DIAMONDTAX      } from '../modules/local/sumtaxonomy/'
include { TIDYVERSE_STRIPCDSPREFIX           } from '../modules/local/tidyverse/stripcdsprefix/'
include { WRITESPADESYAML                    } from '../modules/local/spades/writeyaml/'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { validateInputSamplesheet       } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'
include { typecastBooleanParam           } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'
include { typecastIntegerParam           } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'

//
// SUBWORKFLOW: Consisting of local modules
//
include { EGGNOG                  } from '../subworkflows/local/eggnog/'
include { EUKULELE                } from '../subworkflows/local/eukulele/'
include { HMMCLASSIFY             } from '../subworkflows/local/hmmclassify/'
include { PROKKA_SUBSETS          } from '../subworkflows/local/prokka/subsets/'
include { FASTQC_TRIMGALORE       } from '../subworkflows/local/fastqc/trimgalore/'
include { PRODIGAL                } from '../subworkflows/local/prodigal/'
include { KOFAMSCAN               } from '../subworkflows/local/kofamscan/'
include { DBCAN                   } from '../subworkflows/local/dbcan/'
include { TRANSDECODER            } from '../subworkflows/local/transdecoder/'
include { METAEUK                 } from '../subworkflows/local/metaeuk/'
include { USER_ORFS               } from '../subworkflows/local/user_orfs/'
include { PIPELINE_INITIALISATION } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'
include { PIPELINE_COMPLETION     } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BBMAP_ALIGN                                } from '../modules/nf-core/bbmap/align/'
include { BBMAP_BBDUK                                } from '../modules/nf-core/bbmap/bbduk/'
include { BBMAP_BBNORM                               } from '../modules/nf-core/bbmap/bbnorm/'
include { BBMAP_INDEX                                } from '../modules/nf-core/bbmap/index/'
include { BEDTOOLS_SORT                              } from '../modules/nf-core/bedtools/sort/'
include { SEQKIT_GREP                                } from '../modules/nf-core/seqkit/grep/'
include { CAT_FASTQ            	                     } from '../modules/nf-core/cat/fastq/'
include { CUSTOM_COLLECTFEATURECOUNTS                } from '../modules/nf-core/custom/collectfeaturecounts/main'
include { CUSTOM_COLLECTSTATS                        } from '../modules/nf-core/custom/collectstats/main'
include { DIAMOND_BLASTP as DIAMOND_TAXONOMY         } from '../modules/nf-core/diamond/blastp/'
include { DUCKDB_TABLE2PARQUET                       } from '../modules/nf-core/duckdb/table2parquet/main'
include { FASTQC                                     } from '../modules/nf-core/fastqc/'
include { MEGAHIT                                    } from '../modules/nf-core/megahit/'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/'
include { PIGZ_COMPRESS as PIGZ_ASSEMBLY             } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_DIAMOND_LINEAGE      } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_PE_READS_FWD         } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_PE_READS_REV         } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_SE_READS             } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_BED     } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_CDS     } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_GFF     } from '../modules/nf-core/pigz/compress/'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_PEP     } from '../modules/nf-core/pigz/compress/'
include { PIGZ_UNCOMPRESS as UNPIGZ_GFF              } from '../modules/nf-core/pigz/uncompress/'
include { QUAST                                      } from '../modules/nf-core/quast/'
include { SAMTOOLS_TRIMHEADER                        } from '../modules/nf-core/samtools/trimheader/'
include { SEQTK_MERGEPE                              } from '../modules/nf-core/seqtk/mergepe/'
include { SEQTK_SEQ as SEQTK_SEQ_CONTIG_FILTER       } from '../modules/nf-core/seqtk/seq/'
include { SPADES                                     } from '../modules/nf-core/spades/'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/subread/featurecounts/'
include { TAXONKIT_LINEAGE                           } from '../modules/nf-core/taxonkit/lineage/'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { paramsSummaryMap                           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_SORT_STATS_SAMTOOLS                    } from '../subworkflows/nf-core/bam_sort_stats_samtools/'
include { MMSEQS_FASTA_CLUSTER                       } from '../subworkflows/nf-core/mmseqs_fasta_cluster/'
include { UTILS_NEXTFLOW_PIPELINE                    } from '../subworkflows/nf-core/utils_nextflow_pipeline/'
include { UTILS_NFCORE_PIPELINE                      } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METATDENOVO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_diamond_dbs // channel: paths to Diamond taxonomy databases, read from --diamond_dbs
    ch_user_orfs   // channel: [ meta, gff, faa ], user-provided ORF calls read from --user_orfs
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // CLI params arrive as strings, and 'false' is truthy; nf-schema does not cast them
    def annotate_only_consolidated = typecastBooleanParam('annotate_only_consolidated')
    def save_parquet               = typecastBooleanParam('save_parquet')
    def skip_dbcan                 = typecastBooleanParam('skip_dbcan')
    def skip_eggnog                = typecastBooleanParam('skip_eggnog')
    def skip_eukulele              = typecastBooleanParam('skip_eukulele')
    def skip_fastqc                = typecastBooleanParam('skip_fastqc')
    def skip_kofamscan             = typecastBooleanParam('skip_kofamscan')
    def skip_protein_consolidation = typecastBooleanParam('skip_protein_consolidation')
    def skip_qc                    = typecastBooleanParam('skip_qc')
    def skip_trimming              = typecastBooleanParam('skip_trimming')
    def min_contig_length          = typecastIntegerParam('min_contig_length')
    def trim_bam_header_above      = typecastIntegerParam('trim_bam_header_above')

    if ( ( params.assembler && params.user_assembly ) || ( ! params.assembler && ! params.user_assembly ) ) {
        error "Provide either `--assembler` or `--user_assembly`!"
    }

    if ( params.user_orfs_gff && ! params.user_orfs_faa ) {
        error 'When supplying ORFs via --user_orfs_gff/--user_orfs_faa, both must be specified, --user_orfs_faa file is missing!'
    } else if ( params.user_orfs_faa && ! params.user_orfs_gff ) {
        error 'When supplying ORFs via --user_orfs_gff/--user_orfs_faa, both must be specified, --user_orfs_gff file is missing!'
    }

    // ORF sources are additive, not mutually exclusive
    if ( ! params.orf_caller && ! params.user_orfs && ! ( params.user_orfs_gff && params.user_orfs_faa ) ) {
        error "Provide `--orf_caller`, `--user_orfs`, `--user_orfs_gff`+`--user_orfs_faa`, or a combination of these!"
    }

    // User ORFs refer to contig ids of a fixed assembly, which a fresh --assembler run won't match
    if ( params.assembler && ( params.user_orfs || ( params.user_orfs_gff && params.user_orfs_faa ) ) ) {
        error "You can't input your own ORFs (`--user_orfs`/`--user_orfs_gff`+`--user_orfs_faa`) if you call for assembly with `--assembler`."
    }

    orf_callers = params.orf_caller ? params.orf_caller.tokenize(',').collect { caller -> caller.trim() } : []
    def valid_orf_callers = ['prodigal', 'prokka', 'transdecoder', 'metaeuk']
    orf_callers.each { caller ->
        if ( ! (caller in valid_orf_callers) ) {
            error "Unknown --orf_caller value '${caller}'. Valid values: ${valid_orf_callers.join(', ')}"
        }
    }

    // Fail fast: metaeuk's own error on a non-mmseqs2 directory is easy to mistake for an OOM kill
    if ( 'metaeuk' in orf_callers && params.metaeuk_db ) {
        def metaeuk_db_path = file(params.metaeuk_db)
        if ( metaeuk_db_path.isDirectory() && ! metaeuk_db_path.listFiles().any { f -> f.name.endsWith('.version') } ) {
            error "--metaeuk_db points at a directory (${params.metaeuk_db}) with no '*.version' file inside -- this doesn't look like an mmseqs2-formatted MetaEuk database. See docs/usage.md for the expected layout."
        }
    }

    // Read synchronously so name collisions can `error` before any task runs. User ORF names act as
    // caller names downstream, so must not clash with callers, each other or reserved names.
    def user_orf_names = ( params.user_orfs ? file(params.user_orfs).splitCsv(header: true).collect { row -> row.name } : [] ) +
        ( params.user_orfs_gff && params.user_orfs_faa ? [ params.user_orfs_name ] : [] )
    // Identity in the name keeps runs at different --cluster_min_seq_id apart. BigDecimal, not
    // Math.round, so 0.995 and 1.0 differ: 0.99 -> 99, 1.0 -> 100, 0.995 -> 99_5.
    cluster_pct              = new java.math.BigDecimal(params.cluster_min_seq_id.toString()).multiply(new java.math.BigDecimal("100"))
    protein_consolidate_name = "protein_consolidate_" + cluster_pct.stripTrailingZeros().toPlainString().replace('.', '_')

    def reserved_caller_names = orf_callers + ['locus_consolidate', protein_consolidate_name]
    user_orf_names.each { name ->
        if ( name in reserved_caller_names ) {
            error "--user_orfs/--user_orfs_name '${name}' collides with an active --orf_caller value or a name the pipeline reserves for itself. Pick a different name."
        }
    }
    def duplicate_user_orf_names = user_orf_names.countBy { it }.findAll { _name, count -> count > 1 }.keySet()
    if ( duplicate_user_orf_names ) {
        error "Duplicate user-supplied-ORFs name(s): ${duplicate_user_orf_names.join(', ')}. Every name (--user_orfs rows and --user_orfs_name) must be unique."
    }
    // Caller names become filename components that CUSTOM_COLLECTSTATS splits on dots
    def dotted_user_orf_names = user_orf_names.findAll { name -> name.contains('.') }
    if ( dotted_user_orf_names ) {
        error "--user_orfs/--user_orfs_name name(s) ${dotted_user_orf_names.join(', ')} contain a '.', which is not allowed in a caller name. Pick a different name."
    }

    // Error, not warning: the double-counted cluster table looks ordinary. `toss` is allowed, its
    // counts are conservative, not wrong.
    if ( params.bbmap_ambiguous == 'all' && ! params.featurecounts_fraction && ! skip_protein_consolidation && ( orf_callers || user_orf_names ) ) {
        error "`--bbmap_ambiguous all` counts a multi-mapping read at full weight at every site it aligns to, which double-counts it when protein consolidation sums counts across a cluster. Add `--featurecounts_fraction` so each alignment is weighted 1/N, or `--skip_protein_consolidation` if you do not need the consolidated table."
    }

    assembler     = params.assembler
    assembly_name = params.assembler ?: params.user_assembly_name

    // Display label only, not tied to which callers run
    orfs_name  = params.orf_caller ?: params.user_orfs_name

    ch_hmmrs = channel.empty()
    if ( params.hmmdir ) {
        channel
            .fromPath(params.hmmdir + params.hmmpattern, checkIfExists: true)
            .set { ch_hmmrs }
    } else if ( params.hmmfiles ) {
        channel
            .fromList( params.hmmfiles.tokenize(',') )
            .map { hmmfile -> [ file(hmmfile) ] }
            .set { ch_hmmrs }
    }

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // DL: I'm not sure which parts are still required after nf-schema. The branch { } certainly is needed.
    ch_fastq = ch_samplesheet
        .flatMap { meta, fastq_files ->
            if (fastq_files.size() <= 2) {
                return [[ meta.id, [meta], fastq_files ]]
            } else {
                def pairs = fastq_files.collate(2)
                return [[ meta.id, pairs.collect { pair -> meta + [id: "${meta.id}_${pairs.indexOf(pair) + 1}"] }, fastq_files ]]
            }
        }
        .map { row -> validateInputSamplesheet(row) }
        .branch {
            meta, fastqs ->
                single  : ( meta.single_end && fastqs.size() == 1 ) || ( ! meta.single_end && fastqs.size == 2 )
                    return [ meta, fastqs ]
                multiple: true
                    return [ meta, fastqs ]
        }

    //
    // MODULE: Concatenate FastQ files from the same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    // Only single-row samples: CAT_FASTQ already gzips the rest

    fwd = ch_fastq.single
        .filter { meta, _f -> ! meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[0] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_PE_READS_FWD(fwd.unzipped)

    rev = ch_fastq.single
        .filter { meta, _f -> ! meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[1] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_PE_READS_REV(rev.unzipped)

    se = ch_fastq.single
        .filter { meta, _f -> meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[0] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_SE_READS(se.unzipped)

    ch_fastq = fwd.zipped.concat(PIGZ_PE_READS_FWD.out.archive)
        .join(rev.zipped.concat(PIGZ_PE_READS_REV.out.archive))
        .map { meta, fwd_read, rev_read -> [ meta, [ fwd_read, rev_read ] ] }
        .concat(
            se.zipped
                .concat(PIGZ_SE_READS.out.archive)
                .map { meta, fastq -> [ meta, [ fastq ] ] }
        )
        .concat(CAT_FASTQ.out.reads)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_fastq,
        skip_fastqc || skip_qc,
        skip_trimming
    )

    ch_collect_stats = ch_fastq
        .collect { meta, _fasta -> meta }
        .map { metas -> [ [ id:"${assembly_name}.${orfs_name}" ], metas ] }

    if ( skip_trimming ) {
        ch_collect_stats = ch_collect_stats
            .map { meta, samples -> [ meta, samples, [] ] }

    } else {
        ch_collect_stats = ch_collect_stats
            .combine(
                FASTQC_TRIMGALORE.out.trim_log
                    .collect { _meta, report ->
                        if ( report in List ) {
                            report[0]
                        } else {
                            report
                        }
                    }
                    .map { report -> [ report ] }
            )
    }

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.fastqc_zip.collect { _meta, zip -> zip })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_log.collect { _meta, log -> log })
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect { _meta, zip -> zip })

    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, channel.fromPath(params.sequence_filter).first() )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { _meta, log ->  log }.map { log -> [ log ] }
        ch_collect_stats = ch_collect_stats.combine(ch_bbduk_logs)
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{ _meta, log -> log })
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = channel.empty()
        ch_collect_stats = ch_collect_stats
            .map { meta, samples, report -> [ meta, samples, report, [] ] }
    }

    //
    // MODULE: Interleave sequences for assembly
    //
    ch_interleaved = channel.empty()
    if ( ! params.user_assembly ) {
        SEQTK_MERGEPE(ch_clean_reads)
        ch_interleaved = SEQTK_MERGEPE.out.reads
    }

    //
    // SUBWORKFLOW: Perform digital normalization.
    //
    if ( ! params.user_assembly ) {
        if ( params.bbnorm ) {
            BBMAP_BBNORM(
                ch_interleaved
                    .collect { _meta, fastq -> fastq }
                    .map { fastq -> [ [id:'all_samples', single_end:true], fastq ] }
            )
            ch_pe_reads_to_assembly = BBMAP_BBNORM.out.fastq.map { _meta, fasta -> fasta }
            ch_se_reads_to_assembly = channel.empty()
        } else {
            ch_pe_reads_to_assembly = ch_interleaved
                .filter { meta, _fastq -> ! meta.single_end }
                .map { _meta, fastq -> fastq }
            ch_se_reads_to_assembly = ch_interleaved
                .filter { meta, _fastq -> meta.single_end }
                .map { _meta, fastq -> fastq }
        }
    }

    //
    // MODULE: Run Megahit or Spades on all interleaved fastq files
    //
    if ( params.user_assembly ) {
        // Downstream assumes a gzipped assembly
        if ( ! params.user_assembly.endsWith('.gz') ) {
            PIGZ_ASSEMBLY(
                channel
                    .fromPath(params.user_assembly)
                    .map { path -> [ [ id:params.user_assembly ], path ] }
            )
            ch_assembly_contigs = PIGZ_ASSEMBLY.out.archive.first()
        } else {
            ch_assembly_contigs = channel
                .value ( [ [ id: assembly_name ], file(params.user_assembly) ] )
        }
    } else if ( assembler == 'spades' ) {
        WRITESPADESYAML (
            ch_pe_reads_to_assembly.toList(),
            ch_se_reads_to_assembly.toList()
        )

        ch_spades = ch_pe_reads_to_assembly
            .mix(ch_se_reads_to_assembly)
            .collect()
            .map { it -> [ [ id: assembly_name ], it, [], [] ] }
        SPADES (
            ch_spades,
            WRITESPADESYAML.out.yaml,
            []
        )

        ch_spades_assembly = SPADES.out.transcripts
            .ifEmpty { [] }
            .combine(SPADES.out.contigs.ifEmpty { [] } )

        FORMATSPADES( ch_spades_assembly.first() )
        ch_assembly_contigs = FORMATSPADES.out.assembly
    } else if ( assembler == 'megahit' ) {
        ch_megahit_reads = ch_se_reads_to_assembly.toList()
            .map { se_reads -> [ [ id: 'megahit_assembly', single_end: true ], se_reads, [] ] }

        MEGAHIT(
            ch_megahit_reads,
            ch_pe_reads_to_assembly.toList()
        )
        ch_assembly_contigs = MEGAHIT.out.contigs
            .map { _meta, contigs -> [ [ id: assembly_name ], contigs ] }
    } else {
        error 'Assembler not specified!'
    }

    // Length threshold is set in modules.config
    if ( min_contig_length > 0 ) {
        SEQTK_SEQ_CONTIG_FILTER ( ch_assembly_contigs )
        ch_assembly_contigs = SEQTK_SEQ_CONTIG_FILTER.out.fastx
    }

    //
    // Call ORFs
    //
    ch_gff      = channel.empty()
    ch_protein  = channel.empty()
    // Gathered for a single UNPIGZ_GFF call: an unaliased process can only be invoked once per workflow
    ch_gff_gz   = channel.empty()

    ch_parquet_tables = channel.empty()

    //
    // SUBWORKFLOW: Run Prokka on batches of the assembly
    //
    if ( 'prokka' in orf_callers ) {
        PROKKA_SUBSETS(ch_assembly_contigs, params.prokka_batchsize)
        ch_protein       = ch_protein.mix( PROKKA_SUBSETS.out.faa.map { meta, faa -> [ meta + [caller: 'prokka'], faa ] } )
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA_SUBSETS.out.prokka_log)

        ch_gff_gz = ch_gff_gz.mix( PROKKA_SUBSETS.out.gff.map { meta, gff -> [ meta + [caller: 'prokka'], gff ] } )

        ch_parquet_tables = ch_parquet_tables.mix( PROKKA_SUBSETS.out.gfftsv.map { _meta, tsv -> tsv } )
    }

    //
    // MODULE: Run PRODIGAL on assembly output.
    //
    if ( 'prodigal' in orf_callers ) {
        PRODIGAL( ch_assembly_contigs.map { _meta, contigs -> [ [id: "${assembly_name}.prodigal", caller: 'prodigal'], contigs  ] } )
        ch_protein      = ch_protein.mix(PRODIGAL.out.faa)
        ch_gff_gz       = ch_gff_gz.mix(PRODIGAL.out.gff)
        // No MultiQC module for Prodigal, so write custom content
        ch_multiqc_files = ch_multiqc_files.mix(
            PRODIGAL.out.faa.collectFile { meta, faa ->
                def n_orfs   = 0
                def total_aa = 0
                def reader   = faa.name.endsWith('.gz') ?
                    new java.util.zip.GZIPInputStream(faa.newInputStream()).newReader() :
                    faa.newReader()
                reader.eachLine { line ->
                    if (line.startsWith('>')) {
                        n_orfs = n_orfs + 1
                    } else {
                        total_aa = total_aa + line.trim().length()
                    }
                }
                reader.close()
                def mean_aa = n_orfs > 0 ? (total_aa / n_orfs) : 0.0
                def content = "Sample,n_orfs,total_aa,mean_aa_length\n${meta.id},${n_orfs},${total_aa},${String.format(Locale.ROOT, '%.1f', mean_aa)}\n"
                [ 'prodigal_stats_mqc.csv', content ]
            }
        )
    }

    //
    // SUBWORKFLOW: run TRANSDECODER. Orf caller alternative for eukaryotes.
    //
    if ( 'transdecoder' in orf_callers ) {
        TRANSDECODER (
            ch_assembly_contigs.map { _meta, contigs -> [ [id: "${assembly_name}.transdecoder", caller: 'transdecoder' ], contigs ] },
            params.transdecoder_batchsize
        )
        ch_gff      = ch_gff.mix(TRANSDECODER.out.gff)
        ch_protein  = ch_protein.mix(TRANSDECODER.out.pep)

        PIGZ_TRANSDECODER_BED(TRANSDECODER.out.bed)
        PIGZ_TRANSDECODER_CDS(TRANSDECODER.out.cds)
        PIGZ_TRANSDECODER_GFF(TRANSDECODER.out.gff)
        PIGZ_TRANSDECODER_PEP(TRANSDECODER.out.pep)

        ch_multiqc_files = ch_multiqc_files.mix(
            TRANSDECODER.out.pep.collectFile { meta, pep ->
                def n_orfs   = 0
                def total_aa = 0
                def reader   = pep.name.endsWith('.gz') ?
                    new java.util.zip.GZIPInputStream(pep.newInputStream()).newReader() :
                    pep.newReader()
                reader.eachLine { line ->
                    if (line.startsWith('>')) {
                        n_orfs = n_orfs + 1
                    } else {
                        total_aa = total_aa + line.trim().length()
                    }
                }
                reader.close()
                def mean_aa = n_orfs > 0 ? (total_aa / n_orfs) : 0.0
                def content = "Sample,n_orfs,total_aa,mean_aa_length\n${meta.id},${n_orfs},${total_aa},${String.format(Locale.ROOT, '%.1f', mean_aa)}\n"
                [ 'transdecoder_stats_mqc.csv', content ]
            }
        )
    }

    //
    // SUBWORKFLOW: run METAEUK. Splice-aware ORF caller alternative for eukaryotes.
    //
    if ( 'metaeuk' in orf_callers ) {
        METAEUK (
            ch_assembly_contigs.map { _meta, contigs -> [ [id: "${assembly_name}.metaeuk", caller: 'metaeuk' ], contigs ] },
            params.metaeuk_db,
            params.metaeuk_db_name,
            params.metaeuk_batchsize
        )
        ch_protein = ch_protein.mix(METAEUK.out.faa)

        ch_gff_gz  = ch_gff_gz.mix(METAEUK.out.gff)

        ch_multiqc_files = ch_multiqc_files.mix(
            METAEUK.out.faa.collectFile { meta, faa ->
                def n_orfs   = 0
                def total_aa = 0
                def reader   = faa.name.endsWith('.gz') ?
                    new java.util.zip.GZIPInputStream(faa.newInputStream()).newReader() :
                    faa.newReader()
                reader.eachLine { line ->
                    if (line.startsWith('>')) {
                        n_orfs = n_orfs + 1
                    } else {
                        total_aa = total_aa + line.trim().length()
                    }
                }
                reader.close()
                def mean_aa = n_orfs > 0 ? (total_aa / n_orfs) : 0.0
                def content = "Sample,n_orfs,total_aa,mean_aa_length\n${meta.id},${n_orfs},${total_aa},${String.format(Locale.ROOT, '%.1f', mean_aa)}\n"
                [ 'metaeuk_stats_mqc.csv', content ]
            }
        )
    }

    // User ORF sets act as further callers from here on. ch_gff must stay uncompressed.
    ch_user_orfs_single = params.user_orfs_gff && params.user_orfs_faa ?
        channel.value( [ [ id: params.user_orfs_name ], file(params.user_orfs_gff), file(params.user_orfs_faa) ] ) :
        channel.empty()
    ch_user_orfs_named = ch_user_orfs.mix(ch_user_orfs_single)
        .map { meta, gff, faa -> [ meta + [caller: meta.id, id: "${assembly_name}.${meta.id}"], gff, faa ] }

    USER_ORFS ( ch_user_orfs_named )

    ch_gff_gz = ch_gff_gz.mix( USER_ORFS.out.gff.filter { _meta, gff -> gff =~ /\.gz$/ } )
    ch_gff    = ch_gff.mix( USER_ORFS.out.gff.filter { _meta, gff -> ! (gff =~ /\.gz$/) } )
    ch_protein = ch_protein.mix( USER_ORFS.out.faa )

    UNPIGZ_GFF(ch_gff_gz)
    ch_gff = ch_gff.mix(UNPIGZ_GFF.out.file)

    // Only calls from different callers merge, so one caller gives its own table. Not bedtools merge:
    // overlap alone fuses adjacent genes of one caller.
    FORMAT_GFF2BED ( ch_gff )

    ch_locus_bed = FORMAT_GFF2BED.out.bed
        .map { _meta, bed -> bed }
        .collectFile(name: "${assembly_name}.combined.bed", sort: true)
        .map { bed -> [ [ id: "${assembly_name}.locus_consolidate", caller: 'locus_consolidate' ], bed ] }

    BEDTOOLS_SORT ( ch_locus_bed, [] )
    FORMAT_LOCUSCONSOLIDATE ( BEDTOOLS_SORT.out.sorted )

    ch_gff = ch_gff.mix(FORMAT_LOCUSCONSOLIDATE.out.gff)

    // Cluster locus proteins to join one gene called on different contigs. Members must be loci,
    // since the counts table sums per-locus counts. params.user_orfs is null for a gff/faa pair.
    ch_protein_clusters = channel.empty()
    if ( ! skip_protein_consolidation && ( orf_callers || user_orf_names ) ) {
        // Joined on meta.id to stay 1:1 with several assemblies; sorted for a stable task hash
        FORMAT_LOCUSFAA (
            FORMAT_LOCUSCONSOLIDATE.out.members
                .map { meta, members -> [ meta.id, meta, members ] }
                .join(
                    ch_protein
                        .map { meta, faa -> [ meta.caller, faa ] }
                        .toList()
                        .map { pairs ->
                            def sorted = pairs.sort { a, b -> a[0] <=> b[0] }
                            [ "${assembly_name}.locus_consolidate", sorted.collect { pair -> pair[0] }, sorted.collect { pair -> pair[1] } ]
                        }
                )
                .map { _id, meta, members, callers, faas -> [ meta, members, callers, faas ] }
        )

        MMSEQS_FASTA_CLUSTER (
            FORMAT_LOCUSFAA.out.faa
                .map { _meta, faa -> [ [ id: "${assembly_name}.${protein_consolidate_name}", caller: protein_consolidate_name ], faa ] },
            'linclust'
        )

        // Smallest member as representative, so cluster ids do not depend on MMseqs2 input order
        FORMAT_CLUSTERREPS ( MMSEQS_FASTA_CLUSTER.out.clusters )
        ch_protein_clusters = FORMAT_CLUSTERREPS.out.clusters

        // 'fa', not 'faa': SEQKIT_GREP only declares *.fa.gz/*.fq.gz outputs
        SEQKIT_GREP (
            MMSEQS_FASTA_CLUSTER.out.seqs,
            FORMAT_CLUSTERREPS.out.representatives.map { _meta, representatives -> representatives },
            'fa'
        )

        ch_protein = ch_protein.mix(SEQKIT_GREP.out.filter)
    }

    total_orf_sources = orf_callers.size() + user_orf_names.size()
    if ( annotate_only_consolidated && ! skip_protein_consolidation && total_orf_sources > 1 ) {
        ch_protein = ch_protein.filter { meta, _protein -> meta.caller == protein_consolidate_name }
    }

    //
    // MODULE: Create a BBMap index
    //
    BBMAP_INDEX(ch_assembly_contigs.map { _meta, contigs -> contigs })

    //
    // MODULE: Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, BBMAP_INDEX.out.index )

    //
    // SUBWORKFLOW: classify ORFs with a set of hmm files
    //
    ch_hmmclassify = ch_hmmrs
        .combine(ch_protein)
        .map { hmm, meta, protein -> [ meta, hmm, protein ] }
    HMMCLASSIFY ( ch_hmmclassify )

    //
    // MODULE: FeatureCounts
    //
    BAM_SORT_STATS_SAMTOOLS (
        BBMAP_ALIGN.out.bam,
        ch_assembly_contigs.map { meta, fasta -> [meta, fasta, []] }
    )

    // featureCounts crashes on a BAM header over 2 GiB, i.e. about 75M contigs; idxstats is about as large as the header
    ch_sorted_bam = BAM_SORT_STATS_SAMTOOLS.out.bam
        .join(BAM_SORT_STATS_SAMTOOLS.out.index)
        .join(BAM_SORT_STATS_SAMTOOLS.out.idxstats)
        .branch { _meta, _bam, _bai, idxstats ->
            trim: idxstats.size() > trim_bam_header_above
            keep: true
        }

    SAMTOOLS_TRIMHEADER ( ch_sorted_bam.trim.map { meta, bam, bai, _idxstats -> [ meta, bam, bai ] } )

    ch_featurecounts = SAMTOOLS_TRIMHEADER.out.bam
        .mix( ch_sorted_bam.keep.map { meta, bam, _bai, _idxstats -> [ meta, bam ] } )
        .combine(ch_gff)   // every sample x every caller
        .map { sampleMeta, bam, callerMeta, gff ->
            // Keep sample meta: FEATURECOUNTS_CDS reads single_end from it
            [ sampleMeta + [ id: "${sampleMeta.id}.${callerMeta.caller}", caller: callerMeta.caller ], bam, gff ]
        }

    ch_collect_stats = ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { _meta, idxstats -> idxstats }.map { idxstats -> [ idxstats ] } )

    FEATURECOUNTS_CDS ( ch_featurecounts)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'metatdenovo_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    ch_collect_feature = FEATURECOUNTS_CDS.out.counts
        .map { meta, fc -> [ meta.caller, fc ] }
        .groupTuple()
        .map { caller, fcs -> [ [ id: "${assembly_name}.${caller}", caller: caller ], fcs ] }
        .branch { meta, _fcs ->
            locus_consolidate: meta.caller == 'locus_consolidate'
            other: true
        }

    // Joined on meta.id: .combine() would cross assemblies if several ever run together
    COLLECT_LOCUSCONSOLIDATE (
        ch_collect_feature.locus_consolidate
            .map { meta, fcs -> [ meta.id, meta, fcs ] }
            .join( FORMAT_LOCUSCONSOLIDATE.out.provenance.map { meta, provenance -> [ meta.id, provenance ] } )
            .map { _id, meta, fcs, provenance -> [ meta, fcs, provenance ] }
    )
    ch_versions = ch_versions.mix(COLLECT_LOCUSCONSOLIDATE.out.versions)

    // Keyed on assembly, as the inputs differ in caller. End-anchored strip: string minus removes the
    // first match, which breaks when the assembly name contains the caller name.
    ch_protein_consolidate_counts = channel.empty()
    if ( ! skip_protein_consolidation && ( orf_callers || user_orf_names ) ) {
        COLLECT_PROTEINCONSOLIDATE (
            ch_collect_feature.locus_consolidate
                .map { meta, fcs -> [ meta.id.replaceAll(java.util.regex.Pattern.quote(".${meta.caller}") + '$', ''), fcs ] }
                .join( FORMAT_LOCUSCONSOLIDATE.out.provenance.map { meta, provenance -> [ meta.id.replaceAll(java.util.regex.Pattern.quote(".${meta.caller}") + '$', ''), provenance ] } )
                .join( ch_protein_clusters.map { meta, clusters -> [ meta.id.replaceAll(java.util.regex.Pattern.quote(".${meta.caller}") + '$', ''), clusters ] } )
                .map { assembly, fcs, provenance, clusters ->
                    [ [ id: "${assembly}.${protein_consolidate_name}", caller: protein_consolidate_name ], fcs, provenance, clusters ]
                }
        )
        ch_protein_consolidate_counts = COLLECT_PROTEINCONSOLIDATE.out.counts
    }

    //
    // MODULE: Unassigned_* counts from featureCounts summaries
    //
    ch_collect_summary = FEATURECOUNTS_CDS.out.summary
        .map { meta, summary -> [ meta.caller, summary ] }
        .groupTuple()
        .map { caller, summaries -> [ [ id: "${assembly_name}.${caller}", caller: caller ], summaries ] }
        // locus_consolidate summaries are only used by protein consolidation
        .filter { meta, _summaries ->
            meta.caller != 'locus_consolidate' || ( ! skip_protein_consolidation && ( orf_callers || user_orf_names ) )
        }

    COLLECT_FEATURECOUNTSSUMMARY ( ch_collect_summary )
    ch_versions = ch_versions.mix(COLLECT_FEATURECOUNTSSUMMARY.out.versions)

    // Protein consolidation re-aggregates locus counts, so reuse the locus Unassigned_* counts for it
    ch_unassigned_per_caller = COLLECT_FEATURECOUNTSSUMMARY.out.unassigned
        .branch { meta, _unassigned ->
            locus_consolidate: meta.caller == 'locus_consolidate'
            other: true
        }
    ch_unassigned_protein_consolidate = channel.empty()
    if ( ! skip_protein_consolidation && ( orf_callers || user_orf_names ) ) {
        // Pad to the 4 provenance columns: CUSTOM_COLLECTSTATS needs all fcs files in a call to share columns
        ch_unassigned_protein_consolidate = ch_unassigned_per_caller.locus_consolidate
            .flatMap { _meta, unassigned -> (unassigned instanceof List) ? unassigned : [ unassigned ] }
            .collectFile { f ->
                // Empty (-stub) or header-only files would break lines[1..-1]
                def lines = file(f).readLines()
                def widened = lines.size() < 2
                    ? ( lines ? lines[0] + '\tcallers\tn_calls\tn_loci\tloci\n' : '' )
                    : ( [ lines[0] + '\tcallers\tn_calls\tn_loci\tloci' ] +
                        lines[1..-1].collect { line -> line + '\t\t\t\t' } ).join('\n') + '\n'
                [ file(f).name, widened ]
            }
            .collect()
            .map { widened -> [ [ id: "${assembly_name}.${protein_consolidate_name}", caller: protein_consolidate_name ], widened ] }
    }
    ch_unassigned_per_caller = ch_unassigned_per_caller.other.mix(ch_unassigned_protein_consolidate)

    CUSTOM_COLLECTFEATURECOUNTS ( ch_collect_feature.other )

    // Strips TransDecoder's cds. prefix; a no-op for other callers
    TIDYVERSE_STRIPCDSPREFIX ( CUSTOM_COLLECTFEATURECOUNTS.out.counts )
    ch_versions           = ch_versions.mix(TIDYVERSE_STRIPCDSPREFIX.out.versions)

    // Must hold every annotated caller: CUSTOM_COLLECTSTATS left-joins onto it
    ch_counts_per_caller  = TIDYVERSE_STRIPCDSPREFIX.out.counts.mix(ch_protein_consolidate_counts)
    ch_fcs_for_summary    = ch_counts_per_caller

    ch_merge_tables = channel.empty()

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //
    if ( ! skip_eggnog ) {
        EGGNOG(ch_protein, ch_fcs_for_summary, params.eggnog_db_url, params.eggnog_dmnd_url, params.eggnog_taxa_url)
        ch_merge_tables   = ch_merge_tables.mix ( EGGNOG.out.sumtable )
        ch_parquet_tables = ch_parquet_tables.mix( EGGNOG.out.emappertsv.map { _meta, tsv -> tsv } )
    }

    //
    // SUBWORKFLOW: run kofamscan on the ORF-called amino acid sequences
    //
    if( !skip_kofamscan ) {
        ch_kofamscan = ch_protein.map { meta, protein -> [ meta, protein ] }
        KOFAMSCAN( ch_kofamscan, ch_fcs_for_summary, params.kofam_ko_list_url, params.kofam_profiles_url )
        ch_merge_tables   = ch_merge_tables.mix ( KOFAMSCAN.out.kofamscan_summary )
        ch_parquet_tables = ch_parquet_tables
            .mix( KOFAMSCAN.out.kofam_table_tsv.map { _meta, tsv -> tsv } )
            .mix( KOFAMSCAN.out.kofam_table_uniq.map { _meta, tsv -> tsv } )
    }

    //
    // SUBWORKFLOW: run dbCAN CAZyme annotation on the ORF-called amino acid sequences
    //
    if( !skip_dbcan ) {
        DBCAN( ch_protein, ch_fcs_for_summary )
        ch_merge_tables   = ch_merge_tables.mix ( DBCAN.out.sumtable )
        ch_parquet_tables = ch_parquet_tables.mix( DBCAN.out.cazyme_annotation.map { _meta, tsv -> tsv } )
    }

    //
    // MODULE: QUAST assembly statistics
    //
    QUAST(
        ch_assembly_contigs,
        channel.value([ [:], [] ]),
        channel.value([ [:], [] ])
    )
    ch_multiqc_files = ch_multiqc_files.mix(QUAST.out.results.map { _meta, dir -> dir })

    //
    // SUBWORKFLOW: Eukulele
    //
    if ( ! skip_eukulele ) {
        // file(), not java.io.File, so s3:// paths work
        d = file(params.eukulele_dbpath)
        if ( ! d.exists() ) {
            d.mkdirs()
        }

        // No --eukulele_db means a user-provided database directory
        ch_eukulele_db = channel.empty()
        if ( params.eukulele_db ) {
            ch_eukulele_db = channel
                .of ( params.eukulele_db )
                .map { db -> [ db, file(params.eukulele_dbpath) ] }
        } else {
            ch_eukulele_db = channel.fromPath(params.eukulele_dbpath, checkIfExists: true)
                .map { path -> [ [], path ] }
        }
        ch_eukulele = ch_protein
            .map { meta, protein -> [ [ id: meta.id, caller: meta.caller ], protein ] }
            .combine( ch_eukulele_db )
            .map { meta, fasta, database, directory -> [ [ id: "${meta.id}.${database}", caller: meta.caller ], fasta, database, directory ] }
        EUKULELE(ch_eukulele, ch_fcs_for_summary)

        ch_merge_tables   = ch_merge_tables.mix(EUKULELE.out.taxonomy_summary)
        ch_parquet_tables = ch_parquet_tables.mix( EUKULELE.out.tax.map { _meta, tsv -> tsv } )
    }

    //
    // MODULE: Diamond taxonomy, every caller x every db
    //
    ch_diamond_input = ch_protein.combine( ch_diamond_dbs.map { db -> [ db[0], db[1] ] } )

    DIAMOND_TAXONOMY(
        ch_diamond_input.map { pm, protein, _dm, _db -> [ pm, protein ] },
        ch_diamond_input.map { _pm, _protein, dm, db -> [ dm, db ] },
        102,
        []
    )

    // Not .join(): callers share db names, and join mishandles duplicate keys
    ch_taxonkit_lineage = DIAMOND_TAXONOMY.out.tsv
        .map { it -> [ [ id: "${it[0].id}.${it[0].db}.lineage", db: it[0].db, caller: it[0].caller ], it[1] ] }
        .combine(ch_diamond_dbs)
        .filter { meta, _tsv, dbMeta, _dmnd, _names, _nodes, _ranks, _parse -> meta.db == dbMeta.id }

    TAXONKIT_LINEAGE(
        ch_taxonkit_lineage.map { it -> [ it[0], [], it[1] ] },
        ch_taxonkit_lineage.map { it -> it[4] },
        ch_taxonkit_lineage.map { it -> it[5] }
    )

    PIGZ_DIAMOND_LINEAGE(
        TAXONKIT_LINEAGE.out.tsv
    )

    FORMAT_DIAMOND_TAX_RANKLIST(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .combine(ch_diamond_dbs)
            .filter { archiveMeta, _archiveFile, dbMeta, _dmnd, _names, _nodes, _ranks, _parse -> archiveMeta.db == dbMeta.id }
            .map { archiveMeta, archiveFile, _dbMeta, _dmnd, _names, _nodes, ranks, _parse -> [ [ id: archiveMeta.id - ".lineage" + ".diamond", db: archiveMeta.db, caller: archiveMeta.caller ], archiveFile, ranks ] }
    )

    FORMAT_DIAMOND_TAX_TAXDUMP(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .combine(ch_diamond_dbs.filter { it -> it[5] })
            .filter { archiveMeta, _archiveFile, dbMeta, _dmnd, _names, _nodes, _ranks, _parse -> archiveMeta.db == dbMeta.id }
            .map { archiveMeta, archiveFile, _dbMeta, _dmnd, names, nodes, ranks, _parse -> [ [ id: archiveMeta.id - ".lineage" + ".diamond", db: archiveMeta.db, caller: archiveMeta.caller ], archiveFile, names, nodes, ranks ] }
    )

    // Not .join(): dbs share caller names
    ch_diamondtax_sum_input = FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy
        .combine( ch_fcs_for_summary )
        .filter { meta, _taxonomy, fcsMeta, _fcs -> meta.caller == fcsMeta.caller }
        .map { meta, taxonomy, _fcsMeta, fcs -> [ meta, meta.db, taxonomy, fcs ] }

    SUM_DIAMONDTAX(ch_diamondtax_sum_input, 'diamondtax')

    ch_merge_tables = ch_merge_tables.mix ( SUM_DIAMONDTAX.out.taxonomy_summary )

    //
    // MODULE: Collect statistics from mapping analysis
    //
    // Never call MERGE_TABLES with zero tables: its pivot_wider fails
    MERGE_TABLES (
        ch_merge_tables
            .map { meta, tsv -> [ meta.caller, tsv ] }
            .groupTuple()
            .map { caller, tsvs -> [ [ id: "${assembly_name}.${caller}", caller: caller ], tsvs ] }
    )

    // Left join, so callers without annotation tables still get stats
    ch_fcs_mergetab_per_caller = ch_counts_per_caller
        .map { meta, fcs -> [ meta.caller, meta, fcs ] }
        .join(
            MERGE_TABLES.out.merged_table.map { meta, mergetab -> [ meta.caller, mergetab ] },
            remainder: true
        )
        // Drop right-only entries: a null meta would crash CUSTOM_COLLECTSTATS
        .filter { _caller, meta, _fcs, _mergetab -> meta != null }
        // Optional output: callers without Unassigned_* rows emit nothing
        .join(
            ch_unassigned_per_caller.map { meta, unassigned -> [ meta.caller, unassigned ] },
            remainder: true
        )
        .filter { _caller, meta, _fcs, _mergetab, _unassigned -> meta != null }
        .map { _caller, meta, fcs, mergetab, unassigned ->
            // A single-file glob match collapses to a bare Path rather than a List.
            def unassignedFiles = unassigned == null ? [] : (unassigned instanceof List ? unassigned : [ unassigned ])
            [ meta, fcs, mergetab ?: [], unassignedFiles ]
        }

    ch_collect_stats = ch_collect_stats
        .combine( ch_fcs_mergetab_per_caller )
        .map { _origMeta, samples, trimlogs, bblogs, idxstats, callerMeta, fcs, mergetab, unassigned ->
            // Each fcs file becomes a column, named from its filename between the first two dots
            [ callerMeta, samples, trimlogs, bblogs, idxstats, [ fcs ] + unassigned, mergetab ]
        }

    CUSTOM_COLLECTSTATS(ch_collect_stats)

    //
    // MODULE: Write summary_tables/ as Parquet, keyed on filename since metas collide
    //
    if ( save_parquet ) {
        DUCKDB_TABLE2PARQUET(
            ch_parquet_tables
                .mix( ch_counts_per_caller.map { _meta, tsv -> tsv } )
                .mix( COLLECT_LOCUSCONSOLIDATE.out.counts.map { _meta, tsv -> tsv } )
                .mix( HMMCLASSIFY.out.hmmrank.map { _meta, tsv -> tsv } )
                .mix( FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy.map { _meta, tsv -> tsv } )
                .mix( FORMAT_DIAMOND_TAX_TAXDUMP.out.taxonomy.map { _meta, tsv -> tsv } )
                .mix( CUSTOM_COLLECTSTATS.out.overall_stats.map { _meta, tsv -> tsv } )
                .map { tsv -> [ [ id: tsv.name.replaceAll(/\.tsv(\.gz)?$/, '') ], tsv ] }
        )
    }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'metatdenovo'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
