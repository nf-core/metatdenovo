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
