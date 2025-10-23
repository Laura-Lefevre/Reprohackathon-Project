This repository is a small Snakemake-based pipeline skeleton. The files under `rules/` implement pipeline steps (downloading, trimming, mapping, counting, stat analysis). Many files are placeholders. Use the guidance below to edit, extend, or implement rules and helper scripts.

Key pointers for working in this repo
- High-level architecture: This is a Snakemake pipeline. The top-level `Snakefile` should orchestrate rules from `rules/*.rules`. Rules are currently split by conceptual stage: download (fastq), trimming, mapping, counting and statistical analysis. Several rule files are empty placeholders — fill them with standard Snakemake `rule` blocks when implementing.
- Where to look: `rules/downloading_fastQ.rules` contains an example bash script (uses sra-toolkit `fasterq-dump`). `README.md` notes the `Snakefile` is the entrypoint. `config.yaml` exists but is currently empty — put pipeline parameters (paths, threads, SRR lists, references) there.
- Project conventions
  - Keep each rule file focused on a single logical step (e.g., `trimming.rules` for trimming commands). Follow Snakemake idioms: input, output, params, threads, shell/conda blocks.
  - Prefer calling small helper scripts in `scripts/` for complex steps; add scripts there and reference them from rules.
  - Use absolute or config-driven paths for data directories; avoid hardcoding home-tilde paths in rules (the `downloading_fastQ.rules` script currently uses `~` — prefer `config.yaml` keys).

Developer workflows and commands
- Running pipeline: implement `Snakefile` to include rules and run with snakemake. Typical commands (examples you can surface in rules or README):
  - snakemake --cores N -s Snakefile
  - snakemake --use-conda --cores N -s Snakefile (if using conda environments per rule)
- Installing tools: the example downloader uses `sudo apt install -y sra-toolkit`. Prefer documenting tool dependencies in `containers/` or a reproducible environment (Conda `environment.yaml` or container images under `containers/`).

Patterns and examples from this repo
- Download: `rules/downloading_fastQ.rules` (bash snippet) — it lists SRR identifiers, calls `fasterq-dump --threads $THREADS --progress --outdir "$OUTDIR" "$SRR"`, and cleans temporary directories. When converting this into a Snakemake `rule`:
  - Make SRR identifiers a config-list in `config.yaml`.
  - Use `shell:` with `{threads}` and `{params.outdir}` rather than hard-coded variables.

Edge-cases and expectations for changes
- Many files are empty placeholders. When you implement behavior, also update `README.md` and `config.yaml` so humans and subsequent agents know how to run the pipeline.
- Be explicit about side effects (downloading deletes temporary folders in the example). Prefer safe defaults (do not delete external data unless explicitly requested).

When editing or adding files, follow these small rules for PR-friendly commits
- Add/modify only what is necessary to implement a rule. Keep changes small and self-contained.
- Add a short docstring / top comment in new rule files explaining inputs/outputs and config keys used.

Files to reference while coding
- `Snakefile` (entrypoint — currently empty)
- `rules/downloading_fastQ.rules` (concrete example)
- `rules/*.rules` (placeholders for other stages)
- `scripts/` (place helper scripts here)
- `containers/` (place container images or instructions for reproducible envs)

If something is missing
- If `config.yaml` is empty, add a minimal example keys block (e.g., data_dir, srr_list, threads, reference_fasta) and update `README.md` with run examples.

Ask for clarification when
- The intended target cluster or execution environment is unknown (local machine vs HPC vs container). State the expected environment and any cluster-specific directives.

If you want, I can:
- Convert `rules/downloading_fastQ.rules` into a Snakemake `rule` that reads SRR IDs from `config.yaml` and writes FASTQ outputs.
- Populate `config.yaml` with recommended keys and a small example.

Please review and tell me which section you want expanded or if you want me to convert the downloader to a working Snakemake rule.
