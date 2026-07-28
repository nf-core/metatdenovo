# AGENTS.md

Instructions for AI coding agents working in this repo (nf-core/metatdenovo).

## Before considering any change ready / before opening a PR

Run these, in order:

1. `prek run -a` — pre-commit hooks (prettier, trailing-whitespace, end-of-file-fixer,
   nextflow-lint). Available in the `nf-core` conda env if not on PATH
   (`conda activate nf-core`).
2. `nf-core pipelines lint` (add `--release` when the PR targets `master`) — nf-core
   community pipeline-standards lint. Also in the `nf-core` conda env.
3. `nextflow lint .` — Nextflow "strict syntax" lint. Run it with two Nextflow versions:
   the minimum declared in `nextflow.config` (`nextflowVersion = '!>=X.Y.Z'`) and the
   latest available. Use `NXF_VER=<version> nextflow lint .` to target a version.

## Module conventions

- `modules/nf-core/` holds vendored nf-core modules — don't hand-edit these directly;
  use `nf-core modules update`/`nf-core modules patch` so changes stay tracked as a
  `.diff` file.
- `modules/local/` holds pipeline-specific modules, laid out the same way as vendored
  ones: `main.nf`, `meta.yml`, `environment.yml` per module/subtool directory (e.g.
  `modules/local/seqtk/hmmhitfaas/`).
- Each process declares `conda "${moduleDir}/environment.yml"` and a matching
  Biocontainer/Singularity container, guarded by `task.ext.when`, with a `stub:` block.
- Tool versions are emitted via `eval(...)` into a `topic: versions` channel, not
  collected through a separate `versions.yml`-merging process.
- `subworkflows/local/` groups related modules into subworkflows (e.g. `eggnog`,
  `kofamscan`, `prokka`); mirror that grouping when adding a new tool with multiple
  steps rather than chaining bare modules in the top-level workflow.
