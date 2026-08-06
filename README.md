# CFTK

CFTK is a cfDNA multimodal epigenetic analysis toolkit for processing
cfMethyl-Seq style data and running downstream methylation, fragmentomics,
visualization, modeling, and report workflows. 
For detailed guidance and tutorial, please refer to the [CFTK website](https://chaorongc.github.io/CFTK/index.html)

You can also use the [CFTK model power calculator](https://cftk-model-power.streamlit.app/) before you start to process your cfDNA cohort.


The package is under active development. The `cftk` command uses a compact
project configuration named `cftk_init.json` plus a TSV sample sheet. Existing
legacy nested configurations remain supported.

The model-development power API is included in the Python distribution, but its
large aggregate reference arrays are not. The repository keeps those arrays in
`data/` for the Streamlit app and source-checkout workflows. Installed callers
must provide that directory explicitly or set `CFTK_MODEL_POWER_DATA`.

## Documentation

The documentation website is built with Sphinx and the PyData Sphinx Theme.

Build it locally with:

```bash
python -m pip install -r docs/requirements.txt
python -m sphinx -b html docs docs/_build/html
```

Then open:

```text
docs/_build/html/index.html
```

## Quick Start

Create the pinned core-processing environment and install CFTK from the
checkout:

```bash
mamba env create -f environment.yml
mamba activate cftk
python -m pip install .
```

The default install is intentionally limited to the beginner processing and QC
workflow. Use `.[analysis]` for differential, power, and MESA commands,
`.[fragmentomics]` for WPS/FinaleToolkit workflows, or `.[web]` for the
Streamlit calculator. CFTK reports the required extra when an optional Python
dependency is missing.

Create a project directory and start the guided initializer:

```bash
mkdir example_study
cd example_study
cftk init
```

`cftk init` uses the current directory and discovers only unambiguous
single-lane FASTQ pairs or BAMs. If it creates a `samples.tsv` template, fill in
the explicit `group` and `role` columns and run `cftk init` again. Initialization
validates the selected profile and writes a portable `cftk.lock.json`, then
builds bwa-meth, samtools, and Picard reference companions.

CFTK automatically installs the managed
`twist_human_methylome_hg38` profile version `1.0.0` under
`CFTK_REFERENCE_ROOT` or `~/.cache/cftk/references`. The profile pins an NCBI
GRCh38 no-alt UCSC-ID FASTA and FASTA index, the matching UCSC analysis-set
2bit genome, and the CFTK covered-target BED at an immutable Git commit.

For batch setup, provide the required inputs explicitly:

```bash
cftk init --non-interactive \
  --sample-sheet samples.tsv
```

Reference profiles live under
`<reference-root>/<profile-id>/<version>/manifest.json`. Set
`CFTK_REFERENCE_ROOT` to relocate a project without editing its JSON. See
`examples/` for the compact config, sample sheet, and manifest formats.

The managed downloader stages each
profile, verifies download and installed-file sizes and SHA-256 hashes, validates
profile compatibility, and publishes the version atomically. See the reference
data documentation before using an external registry.

Inspect available commands:

```bash
cftk --help
```

Run the beginner workflow from the initialized project directory:

```bash
cftk run
```

`cftk run` performs a read-only doctor preflight, then runs core processing and
QC with fail-fast stage boundaries. It validates every required artifact and
writes an immutable attempt directory under
`results/provenance/runs/<run-id>/`, including exact commands, tool versions,
expected outputs and figures, stage states, and `run-summary.html`. Preview the
plan without tool probes or computation with `cftk run --dry-run`. Each attempt
also writes `resource-plan.json` and an `evidence/` directory with
stage/artifact/command tables plus sanitized figure previews; `process.cores`
is the total CPU budget and is divided across concurrent sample commands.

Use `cftk doctor` separately for readiness diagnostics or JSON output to CI.
The individual commands remain available for expert workflows:

```bash
cftk doctor --json
cftk --config cftk_init.json process -s 1 2 3 4
```

After a successful core run, use the role-aware downstream planner before
launching differential, fragmentomics, DMR, MESA, and reporting stages:

```bash
cftk plan
cftk analyze --dry-run
cftk analyze
```

The `auto` preset runs descriptive occupancy/WPS/reporting for a one-group
project and adds differential analysis for an explicit two-group control/case
project. Use `cftk plan --preset all` to inspect the complete downstream graph
before `cftk analyze --preset all`. DMR and MESA remain explicit because they
introduce R/external-tool and modeling dependencies. Every attempt records a
preflight report, stage contracts, command mirror, evidence, and HTML summary
under `results/provenance/analysis-runs/<run-id>/`.

The default Twist targeted profile automatically restricts WPS, occupancy, and
DELFI to panel-overlapping reads and regions. Review the recorded scope in the
plan, run manifest, and each scoped stage's
`fragmentomics_scope.json`; use `--fragmentomics-scope genome` only for
validated whole-genome inputs. The sidecar records the target BED identity,
derived interval counts, and the interpretation note.

Step 3 also writes Picard target and alignment metrics. Schema-v2 projects use
the selected profile's covered-target BED. `--target-bed PATH` remains an
expert one-run override; legacy source checkouts retain the bundled fallback.
Picard uses an explicit 8 GB maximum Java heap by default; advanced projects
can tune it with schema-v2 `process.picard_java_memory`.

CFTK appends every external command to
`<output_dir>/results/provenance/commands.jsonl`. Each command is recorded in
full before execution and receives a completion record with its exit status.
This includes reference preparation, raw processing, QC, DMR, and
fragmentomics commands, including commands launched for parallel samples.

Some workflows require external bioinformatics tools and reference files that
are not installed by Python packaging alone. See the documentation for details.

## Roadmap

The fail-safe beginner `cftk run` is implemented. Clean real-data validation
and biological acceptance remain release work tracked in [TODO.md](TODO.md).
Sambamba is the default duplicate-marking tool; advanced users can select
Picard explicitly in the schema-v2 process configuration.

## Model-Power Calculator

Run the calculator from a repository checkout so it can access the aggregate
reference arrays under `data/`:

```bash
python -m pip install ".[web]"
streamlit run apps/model_power_calculator.py
```

The hosted Streamlit deployment currently requires platform authentication, so
this README does not advertise it as a public calculator. Add the public URL
here after the deployment is accessible without authentication.
