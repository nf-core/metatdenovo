/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local
//
include { COLLECT_FEATURECOUNTS              } from '../modules/local/collect/featurecounts/'
include { COLLECT_LOCUS_CONSOLIDATE          } from '../modules/local/collect/locus_consolidate/'
include { COLLECT_PROTEINCONSOLIDATE         } from '../modules/local/collect/proteinconsolidate/'
include { COLLECT_STATS                      } from '../modules/local/collect/stats/'
include { FORMAT_GFF2BED                     } from '../modules/local/format/gff2bed/'
include { FORMAT_LOCUS_CONSOLIDATE           } from '../modules/local/format/locus_consolidate/'
include { FORMAT_LOCUSFAA                    } from '../modules/local/format/locusfaa/'
include { FORMATSPADES                       } from '../modules/local/format/spades/'
include { MERGE_TABLES                       } from '../modules/local/merge/summary/'
include { FORMAT_DIAMOND_TAX_RANKLIST        } from '../modules/local/diamond/format_tax/ranklist/'
include { FORMAT_DIAMOND_TAX_TAXDUMP         } from '../modules/local/diamond/format_tax/taxdump/'
include { SUMTAXONOMY as SUM_DIAMONDTAX      } from '../modules/local/sumtaxonomy/'
include { TRANSRATE                          } from '../modules/nf-core/transrate/'
include { WRITESPADESYAML                    } from '../modules/local/spades/writeyaml/'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { validateInputSamplesheet       } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'

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
include { BEDTOOLS_MERGE                             } from '../modules/nf-core/bedtools/merge/'
include { BEDTOOLS_SORT                              } from '../modules/nf-core/bedtools/sort/'
include { CAT_FASTQ            	                     } from '../modules/nf-core/cat/fastq/'
include { DIAMOND_BLASTP as DIAMOND_TAXONOMY         } from '../modules/nf-core/diamond/blastp/'
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
include { PIGZ_UNCOMPRESS as UNPIGZ_CONTIGS          } from '../modules/nf-core/pigz/uncompress/'
include { PIGZ_UNCOMPRESS as UNPIGZ_GFF              } from '../modules/nf-core/pigz/uncompress/'
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
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    // Exit if the user provides both --assembler and --user_assembly, or both --orf_caller and either of the two --user_orfs,
    // or none of them
    if ( ( params.assembler && params.user_assembly ) || ( ! params.assembler && ! params.user_assembly ) ) {
        error "Provide either `--assembler` or `--user_assembly`!"
    }
    if ( ( params.orf_caller && params.user_orfs_gff ) || ( ! params.orf_caller && ! params.user_orfs_gff ) ) {
        error "Provide either `--orf_caller` or `--user_orfs_gff`/`--user_orfs_faa`!"
    }

    // Exit if the user forgot one of the two `--user_orfs_*`
    if ( params.user_orfs_gff && ! params.user_orfs_faa ) {
        error 'When supplying ORFs, both --user_orfs_gff and --user_orfs_faa must be specified, --user_orfs_faa file is missing!'
    } else if ( params.user_orfs_faa && ! params.user_orfs_gff ) {
        error 'When supplying ORFs, both --user_orfs_gff and --user_orfs_faa must be specified, --user_orfs_gff file is missing!'
    }

    // Exit if the user set params.assembler plus any of params.user_orfs_*
    if ( params.assembler && ( params.user_orfs_gff || params.user_orfs_faa ) ) {
        error "You can't input your own ORFs (`--user_orfs_*`) if you call for assembly with `--assembler`."
    }

    // Split --orf_caller into a list and make sure every entry is a caller we know about
    orf_callers = params.orf_caller ? params.orf_caller.tokenize(',').collect { caller -> caller.trim() } : []
    def valid_orf_callers = ['prodigal', 'prokka', 'transdecoder', 'metaeuk']
    orf_callers.each { caller ->
        if ( ! (caller in valid_orf_callers) ) {
            error "Unknown --orf_caller value '${caller}'. Valid values: ${valid_orf_callers.join(', ')}"
        }
    }

    // Exit if the user asked for metaeuk without also providing a reference database
    if ( 'metaeuk' in orf_callers && ! params.metaeuk_db ) {
        error "When using `--orf_caller metaeuk`, you must also specify `--metaeuk_db`!"
    }

    // Deal with user-supplied assembly to make sure output names are correct
    assembler     = params.assembler
    assembly_name = params.assembler ?: params.user_assembly_name

    // Deal with params from user-supplied ORFs
    orfs_name  = params.orf_caller ?: params.user_orfs_name

    // Name the cross-contig consolidation level after the clustering identity it was produced at,
    // rounded to whole percent, so reruns at a different --cluster_min_seq_id don't overwrite each
    // other's tables -- the same reason assembler and ORF caller are already in output filenames.
    // Computed once here and threaded from here; it's used as a caller name, so it ends up in the
    // counts table's filename and in every annotation table derived from the cluster representatives.
    protein_consolidate_name = "protein_consolidate_${Math.round(params.cluster_min_seq_id * 100)}"


    // If the user supplied hmm files, we will run hmmsearch and then rank the results.
    // Create a channel for hmm files.
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

    // Form a fastq channel from the samplesheet channel
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
        /** DL: In my testing, this fails as entries come in with single meta, multiple read files when appear for single ends
        .map { id, metas, fastq_files ->
            // Ensure single_end is set correctly in meta
            def updatedMetas = metas
                .collect { it + [single_end: (fastq_files.size() / metas.size() == 1)] }
            return [id, updatedMetas, fastq_files]
        }
        **/
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
    //
    // Gzip unzipped read files
    //
    // We're only doing this for samples having a single row in the sample sheet since those with more than one
    // were gzipped by CAT_FASTQ above.
    //

    // Paired end, forward
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

    // Paired end, reverse
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

    // Single end
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

    // Join the three channels with the originally zipped to form a new ch_fastq of the same structure as the original
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
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )

    ch_collect_stats = ch_fastq
        .collect { meta, _fasta -> meta }
        .map { metas -> [ [ id:"${assembly_name}.${orfs_name}" ], metas ] }

    if ( params.skip_trimming ) {
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
        // If the input assembly is not gzipped, do that since all downstream calls assume this
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
        // 1. Write a yaml file for Spades
        WRITESPADESYAML (
            ch_pe_reads_to_assembly.toList(),
            ch_se_reads_to_assembly.toList()
        )

        // 2. Call the module with a channel with all fastq files plus the yaml
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

    // If the user asked for length filtering, perform that with SEQTK_SEQ (the actual length parameter is used in modules.config)
    if ( params.min_contig_length > 0 ) {
        SEQTK_SEQ_CONTIG_FILTER ( ch_assembly_contigs )
        ch_assembly_contigs = SEQTK_SEQ_CONTIG_FILTER.out.fastx
    }

    //
    // Call ORFs
    //
    ch_gff      = channel.empty()
    ch_protein  = channel.empty()
    // PROKKA_SUBSETS/PRODIGAL/METAEUK all emit gzipped GFFs that need UNPIGZ_GFF -- accumulated here
    // and unzipped in ONE call after all three `if` blocks below, rather than one UNPIGZ_GFF call per
    // block: Nextflow doesn't allow the same unaliased process to be invoked more than once in the
    // same workflow context, which broke as soon as two of these callers were active simultaneously
    // (e.g. --orf_caller prokka,prodigal) even though it worked fine for any single caller alone.
    ch_gff_gz   = channel.empty()

    //
    // SUBWORKFLOW: Run PROKKA_SUBSETS on assmebly output, but split the fasta file in chunks of 10 MB, then concatenate and compress output.
    //
    if ( 'prokka' in orf_callers ) {
        PROKKA_SUBSETS(ch_assembly_contigs, params.prokka_batchsize)
        ch_protein       = ch_protein.mix( PROKKA_SUBSETS.out.faa.map { meta, faa -> [ meta + [caller: 'prokka'], faa ] } )
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA_SUBSETS.out.prokka_log)

        ch_gff_gz = ch_gff_gz.mix( PROKKA_SUBSETS.out.gff.map { meta, gff -> [ meta + [caller: 'prokka'], gff ] } )
    }

    //
    // MODULE: Run PRODIGAL on assembly output.
    //
    if ( 'prodigal' in orf_callers ) {
        PRODIGAL( ch_assembly_contigs.map { _meta, contigs -> [ [id: "${assembly_name}.prodigal", caller: 'prodigal'], contigs  ] } )
        ch_protein      = ch_protein.mix(PRODIGAL.out.faa)
        ch_gff_gz       = ch_gff_gz.mix(PRODIGAL.out.gff)
        // No MultiQC-native module for Prodigal (unlike Prokka, handled via
        // PROKKA_SUBSETS.out.prokka_log above), so build a small custom-content CSV of basic
        // ORF/protein stats from its amino-acid FASTA output.
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
                def mean_aa = n_orfs > 0 ? (total_aa / n_orfs) : 0
                def content = "Sample,n_orfs,total_aa,mean_aa_length\n${meta.id},${n_orfs},${total_aa},${String.format(Locale.ROOT, '%.1f', mean_aa)}\n"
                [ 'prodigal_stats_mqc.csv', content ]
            }
        )
    }

    //
    // SUBWORKFLOW: run TRANSDECODER. Orf caller alternative for eukaryotes.
    //
    if ( 'transdecoder' in orf_callers ) {
        TRANSDECODER ( ch_assembly_contigs.map { _meta, contigs -> [ [id: "${assembly_name}.transdecoder", caller: 'transdecoder' ], contigs ] } )
        ch_gff      = ch_gff.mix(TRANSDECODER.out.gff)
        ch_protein  = ch_protein.mix(TRANSDECODER.out.pep)

        PIGZ_TRANSDECODER_BED(TRANSDECODER.out.bed)
        PIGZ_TRANSDECODER_CDS(TRANSDECODER.out.cds)
        PIGZ_TRANSDECODER_GFF(TRANSDECODER.out.gff)
        PIGZ_TRANSDECODER_PEP(TRANSDECODER.out.pep)

        // No MultiQC-native module for TransDecoder either -- same custom-content approach as
        // the Prodigal branch above.
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
                def mean_aa = n_orfs > 0 ? (total_aa / n_orfs) : 0
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
            file(params.metaeuk_db, checkIfExists: true)
        )
        ch_protein = ch_protein.mix(METAEUK.out.faa)

        ch_gff_gz  = ch_gff_gz.mix(METAEUK.out.gff)

        // No MultiQC-native module for MetaEuk either -- same custom-content approach as
        // the Prodigal/TransDecoder branches above.
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
                def mean_aa = n_orfs > 0 ? (total_aa / n_orfs) : 0
                def content = "Sample,n_orfs,total_aa,mean_aa_length\n${meta.id},${n_orfs},${total_aa},${String.format(Locale.ROOT, '%.1f', mean_aa)}\n"
                [ 'metaeuk_stats_mqc.csv', content ]
            }
        )
    }

    // Single UNPIGZ_GFF call covering every gzipped-GFF caller active this run -- see the ch_gff_gz
    // comment above for why this can't be one call per caller branch.
    UNPIGZ_GFF(ch_gff_gz)
    ch_gff = ch_gff.mix(UNPIGZ_GFF.out.file)

    //
    // Consolidate overlapping same-contig CDS calls from different active callers into single loci
    // before counting, so a read supporting one real gene isn't counted once per caller that called
    // it (#463). Runs unconditionally, even with a single active caller: with nothing to overlap,
    // every bedtools-merge row has exactly one contributor, so FORMAT_LOCUS_CONSOLIDATE's
    // ID-inheritance rule reproduces that caller's own GFF byte-for-byte -- graceful degradation,
    // no separate code path. Rides through the existing, unmodified ch_featurecounts cross-product
    // and FEATURECOUNTS_CDS call below by being mixed into ch_gff as just another caller.
    //
    // versions_gzip/versions_bedtools are topic:versions tuple emits, not versions.yml Paths --
    // like every other module's topic-based version emit in this pipeline, they're picked up
    // automatically by the global channel.topic("versions") collection below and must not be
    // mixed into ch_versions directly (that expects Path entries only).
    FORMAT_GFF2BED ( ch_gff )

    ch_locus_bed = FORMAT_GFF2BED.out.bed
        .map { _meta, bed -> bed }
        .collectFile(name: "${assembly_name}.combined.bed", sort: true)
        .map { bed -> [ [ id: "${assembly_name}.locus_consolidate", caller: 'locus_consolidate' ], bed ] }

    BEDTOOLS_SORT ( ch_locus_bed, [] )
    BEDTOOLS_MERGE ( BEDTOOLS_SORT.out.sorted )
    FORMAT_LOCUS_CONSOLIDATE ( BEDTOOLS_MERGE.out.bed )

    ch_gff = ch_gff.mix(FORMAT_LOCUS_CONSOLIDATE.out.gff)

    // Populate channels if the user provided the orfs
    if ( params.user_orfs_faa && params.user_orfs_gff ) {
        ch_gff = channel.value ( [ [ id: "${assembly_name}.${orfs_name}", caller: orfs_name ], file(params.user_orfs_gff) ] )
        ch_protein = channel.value ( [ [ id: "${assembly_name}.${orfs_name}", caller: orfs_name ], file(params.user_orfs_faa) ] )
    }

    //
    // Consolidate calls for the same gene that ended up on DIFFERENT contigs, which coordinates
    // cannot detect: a splice-aware genomic call and a transcript-derived call for one gene share no
    // coordinate system, but converge on nearly the same protein (#460). Cluster the proteins and
    // treat one cluster as one gene.
    //
    // Clustering runs on loci rather than on each caller's raw ORFs, because the counts table this
    // feeds aggregates the per-locus counts locus consolidation produces above -- so cluster members
    // have to BE loci. FORMAT_LOCUSFAA resolves each locus back to one protein sequence for that.
    //
    // Deliberately placed before ch_protein is consumed below, and dependent only on the GFFs and
    // proteins rather than on any counts, so it runs alongside mapping instead of behind it.
    // Skipped entirely for user-supplied ORFs, which replace ch_gff/ch_protein wholesale just above
    // and already bypass locus consolidation.
    ch_protein_clusters = channel.empty()
    if ( ! params.skip_protein_consolidation && orf_callers ) {
        // Keyed on the locus-consolidation meta.id rather than .combine()d, so this stays a genuine
        // 1:1 pairing if the pipeline ever consolidates more than one assembly in a run. Callers are
        // sorted so the two lists stay aligned and the task's inputs hash reproducibly.
        FORMAT_LOCUSFAA (
            FORMAT_LOCUS_CONSOLIDATE.out.members
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

        ch_protein_clusters = MMSEQS_FASTA_CLUSTER.out.clusters
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
    // MODULE: FeatureCounts. Create a table for each samples that provides raw counts as result of the alignment.
    //
    BAM_SORT_STATS_SAMTOOLS (
        BBMAP_ALIGN.out.bam,
        ch_assembly_contigs.map { meta, fasta -> [meta, fasta, []] }
    )

    ch_featurecounts = BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(ch_gff)   // deliberate cross product: every sample x every active caller
        .map { sampleMeta, bam, callerMeta, gff ->
            // sampleMeta + [...] (not a fresh map) so fields FEATURECOUNTS_CDS itself reads off the
            // sample side (e.g. single_end, which picks the -p/paired-end flag) survive the rename.
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

    // Locus-consolidated counts get their own collect module (provenance join), branched out here
    // so COLLECT_FEATURECOUNTS itself -- and every other caller's table -- stays untouched.
    // Keyed by meta.id (not .combine(), which would cross-join every assembly's featureCounts
    // group against every assembly's provenance file if this pipeline ever ran multiple
    // simultaneous assemblies) -- matches ch_featurecounts/ch_collect_feature's own convention
    // just above of staying joinable rather than relying on positional/cardinality pairing.
    COLLECT_LOCUS_CONSOLIDATE (
        ch_collect_feature.locus_consolidate
            .map { meta, fcs -> [ meta.id, meta, fcs ] }
            .join( FORMAT_LOCUS_CONSOLIDATE.out.provenance.map { meta, provenance -> [ meta.id, provenance ] } )
            .map { _id, meta, fcs, provenance -> [ meta, fcs, provenance ] }
    )
    ch_versions = ch_versions.mix(COLLECT_LOCUS_CONSOLIDATE.out.versions)

    // Third consolidation level: sum the per-locus counts above across each protein cluster, so a
    // gene called on two contigs is reported once. Safe to sum after the fact rather than recount,
    // because a read only aligns to one contig -- provided each read was exclusively assigned to one
    // feature at counting time, which is what --bbmap_ambiguous/--featurecounts_fraction control
    // (#464).
    //
    // Keyed on the assembly name rather than meta.id, since the three inputs carry three different
    // caller names; meta.id is "<assembly>.<caller>" throughout, so stripping the caller recovers a
    // key they share without assuming there is only ever one assembly in flight.
    ch_protein_consolidate_counts = channel.empty()
    if ( ! params.skip_protein_consolidation && orf_callers ) {
        COLLECT_PROTEINCONSOLIDATE (
            ch_collect_feature.locus_consolidate
                .map { meta, fcs -> [ meta.id - ".${meta.caller}", fcs ] }
                .join( FORMAT_LOCUS_CONSOLIDATE.out.provenance.map { meta, provenance -> [ meta.id - ".${meta.caller}", provenance ] } )
                .join( ch_protein_clusters.map { meta, clusters -> [ meta.id - ".${meta.caller}", clusters ] } )
                .map { assembly, fcs, provenance, clusters ->
                    [ [ id: "${assembly}.${protein_consolidate_name}", caller: protein_consolidate_name ], fcs, provenance, clusters ]
                }
        )
        ch_protein_consolidate_counts = COLLECT_PROTEINCONSOLIDATE.out.counts
    }

    COLLECT_FEATURECOUNTS ( ch_collect_feature.other )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    // [ meta(caller), tsv ] -- kept as a proper tuple (not stripped) so every consumer below can
    // join it against other per-caller channels instead of relying on positional channel pairing.
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts

    // Initialize ch_merge_tables that will be populated with tables from annotation tools and used by the MERGE_TABLES module which output will then be passed to the COLLECT_STATS module
    ch_merge_tables = channel.empty()

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //
    if ( ! params.skip_eggnog ) {
        EGGNOG(ch_protein, ch_fcs_for_summary)
        ch_merge_tables = ch_merge_tables.mix ( EGGNOG.out.sumtable )
    }

    //
    // SUBWORKFLOW: run kofamscan on the ORF-called amino acid sequences
    //
    if( !params.skip_kofamscan ) {
        ch_kofamscan = ch_protein.map { meta, protein -> [ meta, protein ] }
        KOFAMSCAN( ch_kofamscan, ch_fcs_for_summary, params.kofam_ko_list_url, params.kofam_profiles_url )
        ch_merge_tables = ch_merge_tables.mix ( KOFAMSCAN.out.kofamscan_summary )
    }

    //
    // SUBWORKFLOW: run dbCAN CAZyme annotation on the ORF-called amino acid sequences
    //
    if( !params.skip_dbcan ) {
        DBCAN( ch_protein, ch_fcs_for_summary )
        ch_merge_tables = ch_merge_tables.mix ( DBCAN.out.sumtable )
    }

    // set up contig channel to use in TransRate
    UNPIGZ_CONTIGS(ch_assembly_contigs)
    ch_unzipped_contigs = UNPIGZ_CONTIGS.out.file

    //
    // MODULE: Use TransRate to judge assembly quality, piped into MultiQC
    //
    TRANSRATE(ch_unzipped_contigs.map { meta, contigs -> [ meta, contigs, [] ] })
    ch_multiqc_files = ch_multiqc_files.mix(
        TRANSRATE.out.assemblies.collectFile { meta, csv -> [ "${meta.id}_assemblies_mqc.csv", csv.text ] }
    )

    //
    // SUBWORKFLOW: Eukulele
    //
    if ( ! params.skip_eukulele ) {
        // Make sure the eukulele_dbpath exists
        d = new File("${params.eukulele_dbpath}")
        if ( ! d.exists() ) {
            d.mkdirs()
        }

        // Create a channel for EUKulele either with a named database or not. The latter means a user-provided database in a directory.
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

        ch_merge_tables = ch_merge_tables.mix(EUKULELE.out.taxonomy_summary)
    }

    //
    // Call Diamond for taxonomy with amino acid sequences
    //
    // Explicit .combine() (rather than relying on ch_protein having historically carried only one
    // item) so every active caller's proteins are searched against every diamond db -- with N>1
    // callers, ch_protein is a genuine multi-item queue channel and an implicit pairing would risk
    // silently degrading to positional lockstep instead of the intended full cross product.
    ch_diamond_input = ch_protein.combine( ch_diamond_dbs.map { db -> [ db[0], db[1] ] } )

    DIAMOND_TAXONOMY(
        ch_diamond_input.map { pm, protein, _dm, _db -> [ pm, protein ] },
        ch_diamond_input.map { _pm, _protein, dm, db -> [ dm, db ] },
        102,
        []
    )

    // Create a unified channel of the output from Diamond together with the diamond db info to
    // make sure the channels are synchronized before calling TAXONKIT_LINEAGE.
    // .join() is unsafe here: with N>1 callers, multiple items (one per caller) share the same db
    // name, but ch_diamond_dbs has only one row per db -- Nextflow's join() silently drops all but
    // one matching left-side item per key when the left side has duplicate keys. .combine() + a
    // .filter() on matching db name gives the same pairing without that pitfall.
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

    // Same .join()-with-duplicate-keys pitfall as above: FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy has
    // multiple items (one per db) sharing the same caller, but ch_fcs_for_summary has one item per
    // caller -- .combine() + .filter() instead of .join().
    ch_diamondtax_sum_input = FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy
        .combine( ch_fcs_for_summary )
        .filter { meta, _taxonomy, fcsMeta, _fcs -> meta.caller == fcsMeta.caller }
        .map { meta, taxonomy, _fcsMeta, fcs -> [ meta, meta.db, taxonomy, fcs ] }

    SUM_DIAMONDTAX(ch_diamondtax_sum_input, 'diamondtax')

    ch_merge_tables = ch_merge_tables.mix ( SUM_DIAMONDTAX.out.taxonomy_summary )

    //
    // MODULE: Collect statistics from mapping analysis
    //
    // MERGE_TABLES only runs for callers that actually have >=1 annotation table -- matching the
    // pre-multi-caller behavior, where MERGE_TABLES was never invoked at all (not invoked-with-an-
    // empty-list; genuinely zero invocations, since .collect() on a channel that never emits produces
    // no output either) whenever every annotation subworkflow was skipped. MERGE_TABLES's own R script
    // isn't written to handle being called with zero input tables (glob-of-nothing breaks its
    // pivot_wider), so preserving "don't call it at all" here, rather than "call it with []", matters.
    MERGE_TABLES (
        ch_merge_tables
            .map { meta, tsv -> [ meta.caller, tsv ] }
            .groupTuple()
            .map { caller, tsvs -> [ [ id: "${assembly_name}.${caller}", caller: caller ], tsvs ] }
    )

    // COLLECT_STATS must still run once per ACTIVE caller regardless of whether that caller got a
    // MERGE_TABLES output -- left-join (remainder: true) onto COLLECT_FEATURECOUNTS.out.counts (which
    // always has exactly one item per active caller) so a caller with no annotation tables still gets
    // a COLLECT_STATS invocation, with mergetab defaulting to [] -- COLLECT_STATS's own script already
    // handles a missing mergetab gracefully (`if (mergetab) {...} else {...}`), mirroring what the
    // pre-multi-caller code's `.ifEmpty { [ [] ] }` fallback did for the single-run case.
    ch_fcs_mergetab_per_caller = COLLECT_FEATURECOUNTS.out.counts
        .map { meta, fcs -> [ meta.caller, meta, fcs ] }
        .join(
            MERGE_TABLES.out.merged_table.map { meta, mergetab -> [ meta.caller, mergetab ] },
            remainder: true
        )
        .map { _caller, meta, fcs, mergetab -> [ meta, fcs, mergetab ?: [] ] }

    ch_collect_stats = ch_collect_stats
        .combine( ch_fcs_mergetab_per_caller )
        .map { _origMeta, samples, trimlogs, bblogs, idxstats, callerMeta, fcs, mergetab ->
            [ callerMeta, samples, trimlogs, bblogs, idxstats, fcs, mergetab ]
        }

    COLLECT_STATS(ch_collect_stats)

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
