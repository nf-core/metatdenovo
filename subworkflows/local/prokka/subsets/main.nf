//
// Run PROKKA on contigs that are split by size, then concatenate output and gunzip it
//

include { PROKKA                      } from '../../../../modules/nf-core/prokka/main'
include { FIND_CONCATENATE as GFF_CAT } from '../../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as FAA_CAT } from '../../../../modules/nf-core/find/concatenate/main'
include { FIND_CONCATENATE as FFN_CAT } from '../../../../modules/nf-core/find/concatenate/main'
include { PROKKAGFF2TSV               } from '../../../../modules/local/prokka/gff2tsv/main'

workflow PROKKA_SUBSETS {
    take:
    contigs   // channel:  tuple val(meta), file(contigs)
    batchsize // channel: strings like '10.MB'. Usually from params.prokka_batchsize

    main:

    PROKKA(
        contigs
            .map { _meta, ctg -> ctg }
            .splitFasta(size: batchsize, file: true)
            .map { ctg -> [ [ id: ctg.getBaseName() ], ctg] },
        [],
        []
    )

    // PROKKA runs once per contig batch (see splitFasta above), so on any assembly larger
    // than params.prokka_batchsize there are multiple *.txt summaries, each with the same
    // literal "organism: Genus species strain " line (no --genus/--species/--strain is
    // passed to PROKKA). MultiQC's prokka module derives the sample name from that line, so
    // feeding the raw per-batch files in directly makes every batch collide under the same
    // "strain" sample name. Sum the numeric fields across all batches into a single
    // per-assembly summary instead, with a real sample name.
    ch_log = contigs
        .map { meta, _contigs -> meta.id }
        .combine(PROKKA.out.txt.map { _meta, txt -> txt }.collect().map { txts -> [ txts ] })
        .map { assembly_id, txts ->
            def totals = [:]
            def order  = []
            txts.each { txt ->
                txt.readLines().each { line ->
                    def parts = line.split(':', 2)
                    if (parts.size() == 2) {
                        def key   = parts[0].trim()
                        def value = parts[1].trim()
                        if (key != 'organism' && value.isInteger()) {
                            if (!totals.containsKey(key)) { order << key }
                            totals[key] = (totals[key] ?: 0) + value.toInteger()
                        }
                    }
                }
            }
            def content = "organism: Genus species ${assembly_id}_prokka\n" +
                order.collect { key -> "${key}: ${totals[key]}" }.join('\n') + '\n'
            // No "_mqc" in this filename: that suffix routes a file into MultiQC's separate
            // custom-content mechanism instead of the normal per-module log search that the
            // native prokka module uses (which detects files by content, not filename).
            [ "${assembly_id}.prokka_summary.txt", content ]
        }
        .collectFile()

    ch_gff = contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.gff.collect { _meta, gff -> gff }.map { gff -> [ gff ] })

    GFF_CAT(ch_gff)

    ch_faa = contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.faa.collect { _meta, protein -> protein }.map { protein -> [ protein ] })

    FAA_CAT(ch_faa)

    ch_ffn = contigs.map{ meta, _contigs -> [ id:"${meta.id}.prokka" ] }
        .combine(PROKKA.out.ffn.collect { _meta, fnn -> fnn }.map { fnn -> [ fnn ] })

    FFN_CAT(ch_ffn)

    PROKKAGFF2TSV(GFF_CAT.out.file_out)

    emit:
    gff        = GFF_CAT.out.file_out.first()
    faa        = FAA_CAT.out.file_out.first()
    ffn        = FFN_CAT.out.file_out.first()
    gfftsv     = PROKKAGFF2TSV.out.tsv
    prokka_log = ch_log
}
