# nf-core/metatdenovo: Coping with large datasets

## Introduction

Large projects -- many samples, deep sequencing, or both -- can push Megahit past the memory available on a given system.
This page gives concrete starting points for three params that can bring memory use back down, in the order we'd suggest trying them, with the caveat below.

Before any of that: the choice of assembler itself matters most.
Megahit (the pipeline's default, `--assembler megahit`) is substantially less memory-hungry than Spades (`--assembler spades`) for the same data, so if you're running Spades and hitting memory limits, switching to Megahit is the single biggest lever available -- try that before reaching for any of the three params below.

> [!WARNING]
> Only one of the numbers below (`--bbnorm_target 50`) comes from real, reported experience.
> The rest are either Megahit's own documented presets or educated-guess starting points, not yet backed by our own benchmarks.
> Treat them as a first attempt, not a tuned recommendation, and please report back what worked (or didn't) for your dataset via a [GitHub issue](https://github.com/nf-core/metatdenovo/issues) -- that's how this page will get more precise over time.

## The order we'd suggest trying these in

1. **Digital normalization** (`--bbnorm`) -- reduces the total volume of read data going into the assembler, so it's the biggest single lever and the one to reach for first if you're far over budget.
2. **`--megahit_min_count`** -- a smaller, more surgical adjustment to Megahit's own graph construction; try this if normalization alone doesn't get you far enough, or if you'd rather not touch the input reads at all.
3. **`--megahit_k_min` / `--megahit_k_max` / `--megahit_k_step`** -- the most aggressive option, since it skips Megahit's most sensitive (and most memory-hungry) low-k iterations entirely. Treat this as a last resort.

See [Usage: Digital normalization](usage.md#digital-normalization) and [Usage: Assembler options](usage.md#assembler-options) for what each param does; this page is just about what values to try.

### 1. Digital normalization

Defaults are `--bbnorm_target 100` and `--bbnorm_min 5`.
In one >90-sample project combining metagenomics and metatranscriptomics, the default target of 100 was not enough to bring the assembly within reach -- lowering `--bbnorm_target` to `50` was needed to get close to a passing assembly.

Suggested ladder: try the default (100) first; if that's not enough, try 50; if still not enough, go lower still, expecting progressively more loss of low-abundance signal as you do.
We don't yet have a confirmed value below 50.
`--bbnorm_min` hasn't been tuned in practice yet -- lowering it (below the default of 5) would retain more of the low-abundance tail rather than discarding it, which is the opposite direction from lowering `--bbnorm_target`, so the two params can be adjusted somewhat independently.

### 2. `--megahit_min_count`

Megahit's own default is 2.
Megahit's presets don't offer a "large dataset" value for this one -- its `meta-sensitive` preset actually lowers it to 1 (for more sensitivity on small assemblies, the opposite goal).
As a starting hint, try 3 as a first, modest step up; we don't have a confirmed value for large datasets yet.

### 3. `--megahit_k_min` / `--megahit_k_max` / `--megahit_k_step`

Megahit itself documents a `meta-large` preset for this scenario -- `--k-min 27 --k-max 127 --k-step 10` -- described in its own `--help` text as intended for "large & complex metagenomes, like soil".
That's the best starting point we have, since it comes from Megahit's own authors rather than from our own testing.
Expect a real trade-off in sensitivity to low-coverage or short reads at these settings -- this is why it's the last lever to reach for, not the first.

## Recovering a Megahit run that was killed partway through

A large assembly can run for days, and Megahit steps up through an increasing series of k-mer sizes as it goes -- so a walltime limit, an out-of-memory kill, a spot-instance eviction, or any other hard stop can land partway through, after real progress has already been made.

`-resume` on its own does not help here. Nextflow's resume works at the level of whole tasks: a task is either cached (it previously finished with exit code 0) or it is not, and Nextflow has no visibility into how far a killed task's own process got internally. A killed `MEGAHIT` task is simply not cached, so a plain `-resume` reruns it from scratch -- discarding however many days of progress it had already made, even though nothing downstream needs to be redone.

Megahit itself, independent of Nextflow, keeps its own checkpoint inside its output directory and can pick up from the last one it completed with `megahit --continue`. That's the tool this section uses to recover, before handing the finished assembly back to the pipeline as a user-provided assembly (see [Assembler options](usage.md#assembler-options)).

> [!WARNING]
> Don't run `nextflow clean` on the run, and don't move or delete the work directory, before completing the recovery below.
> Megahit records every input read file's absolute path at the time it started, and expects the exact same paths (the ones Nextflow staged as symlinks in the task's own work directory) to still resolve when it resumes.

### 1. Find the failed task's work directory

If the run just failed, Nextflow prints it directly in the error output:

```console
Error executing process > 'NFCORE_METATDENOVO:METATDENOVO:MEGAHIT (megahit_assembly)'
...
Work dir:
  /path/to/work/3f/8a1b2c...
```

If you're coming back to this later, or the run is still going and you want to check without stopping it, look it up instead:

```bash
nextflow log <run-name-or-session-id> -f hash,process,status,workdir | grep MEGAHIT
```

(`nextflow log` with no arguments lists recent run names/session IDs.)

### 2. Confirm there's a checkpoint to resume from

```bash
cd /path/to/work/3f/8a1b2c...
ls megahit_out/
```

When running, nf-core/metatdenovo doesn't set Megahit's `-o` explicitly, so Megahit always writes to `megahit_out/` inside the task's own work directory. If you see `checkpoints.txt`, `options.json` and a populated `intermediate_contigs/`, there's a checkpoint to resume from. `tail megahit_out/log` shows how far it got -- e.g. a last line like `Extracting solid (k+1)-mers and building sdbg for k = 69` means every earlier k was completed and only that in-progress step needs to be redone.

If `checkpoints.txt` doesn't exist yet, Megahit hadn't reached its first checkpoint -- there's nothing to resume, and a fresh run is the only option (in practice rare, since the first checkpoint comes early relative to the k-mer stages that take days).

### 3. Make sure the recovery environment has enough resources

`megahit --continue` re-reads every original setting -- including `-t`/`-m` -- from `megahit_out/options.json` and ignores anything passed alongside `--continue` on the command line, so you can't scale threads/memory down for the resume. Whatever machine you run it on needs at least as much CPU and RAM as the original task had (the pipeline's `process_high` label, `conf/base.config` -- 12 CPUs / 72 GB by default, more if the task had already retried).

### 4. Resume Megahit itself, from inside the work directory

Use the exact same container the task used, so the Megahit version matches -- find it in that same work directory's `.command.run` (search for the image reference, e.g. `grep -o '"[^"]*megahit[^"]*"' .command.run`, or for whichever container line your executor prints).

Docker:

```bash
docker run --rm -v "$(pwd)":/data -w /data <image-from-.command.run> megahit --continue -o megahit_out
```

Singularity/Apptainer:

```bash
singularity exec -B "$(pwd)":/data --pwd /data <image-from-.command.run> megahit --continue -o megahit_out
```

Or, if you have the exact same Megahit version available natively (e.g. via conda), just run `megahit --continue -o megahit_out` directly -- no container needed.

Megahit logs `passing check point N` for every step it skips before it gets back to real work, so you'll see immediately whether it picked up where it left off rather than starting over.

### 5. Hand the finished assembly back to the pipeline

Once it finishes, the assembly is `megahit_out/<prefix>.contigs.fa` (uncompressed -- the pipeline's own `.gz` step hasn't run yet, and doesn't need to; `--user_assembly` accepts either). Copy it somewhere durable, then start a **new** pipeline run pointing at it instead of `--assembler megahit`:

```bash
nextflow run nf-core/metatdenovo -profile docker --outdir results/ --input samples.csv \
    --user_assembly /path/to/megahit_assembly.contigs.fa \
    --user_assembly_name megahit_assembly
```

This is a fresh run, not a Nextflow-level `-resume` of the original session -- and that's deliberate. It's possible in principle to hand-place a `.exitcode` and the expected output files into the original task's work directory and `-resume` the same session, but that relies on Nextflow's cache bookkeeping matching by hand, and a small mistake there is easy to get wrong and hard to notice until much later. `--user_assembly` is the pipeline's own tested entry point for supplying a pre-built assembly, and it skips reads-to-assembly steps (`SEQTK_MERGEPE`, digital normalization, the assembler itself) automatically, so nothing upstream of the assembly is redone.
