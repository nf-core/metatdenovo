//
// Normalise user-supplied ORFs (--user_orfs_gff/--user_orfs_faa) before they join the pipeline's
// own callers. A GFF/FASTA pair re-supplied from a MetaEuk run (e.g. re-fed after MetaEuk's own
// resume/cache broke on a large run) is raw MetaEuk output, not yet through the same
// FORMAT_METAEUK_GFF/FORMAT_METAEUKFAA rewrite a pipeline-internal `--orf_caller metaeuk` run
// always gets first -- so it carries `TCS_ID=`/`Target_ID=` attributes instead of a plain `ID=`.
//

include { FORMAT_METAEUKFAA  } from '../../../modules/local/format/metaeukfaa/main'
include { FORMAT_METAEUK_GFF } from '../../../modules/local/format/metaeuk/main'

// A plain ID= attribute is unique to already-normalised GFFs; MetaEuk's own raw output never has
// one (only Target_ID=/TCS_ID=). Reads at most 20 lines, decompressing if needed, regardless of
// how large the underlying file is -- these can be tens of GB for a real assembly.
def isRawMetaeukGff(gff) {
    def stream = gff.name.endsWith('.gz') ?
        new java.util.zip.GZIPInputStream(gff.newInputStream()) :
        gff.newInputStream()
    def reader = new BufferedReader(new InputStreamReader(stream))
    reader.withCloseable { r ->
        r.lines().limit(20).anyMatch { line -> line.contains('TCS_ID=') }
    }
}

workflow USER_ORFS {

    take:
    user_orfs // channel: [ val(meta), path(gff), path(faa) ]

    main:

    ch_branched = user_orfs
        .branch { _meta, gff, _faa ->
            metaeuk_shaped: isRawMetaeukGff(gff)
            generic: true
        }

    FORMAT_METAEUK_GFF ( ch_branched.metaeuk_shaped.map { meta, gff, _faa -> [ meta, gff ] } )
    FORMAT_METAEUKFAA  ( ch_branched.metaeuk_shaped.map { meta, _gff, faa -> [ meta, faa ] } )

    ch_formatted = FORMAT_METAEUK_GFF.out.format_gff
        .join( FORMAT_METAEUKFAA.out.format_faa )
        .mix( ch_branched.generic )

    emit:
    gff = ch_formatted.map { meta, gff, _faa -> [ meta, gff ] } // channel: [ val(meta), path(gff) ]
    faa = ch_formatted.map { meta, _gff, faa -> [ meta, faa ] } // channel: [ val(meta), path(faa) ]
}
