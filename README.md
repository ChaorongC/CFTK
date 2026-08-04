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

Install the Python package from the checkout:

```bash
python -m pip install .
```

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

Run raw processing after initialization:

```bash
cftk --config cftk_init.json process -s 1 2 3 4
```

Step 3 also writes Picard target and alignment metrics. Schema-v2 projects use
the selected profile's covered-target BED. `--target-bed PATH` remains an
expert one-run override; legacy source checkouts retain the bundled fallback.

Some workflows require external bioinformatics tools and reference files that
are not installed by Python packaging alone. See the documentation for details.

## Roadmap

`cftk doctor`, a fail-safe beginner `cftk run`, and real-data end-to-end
validation remain explicit release TODOs. Their prerequisites and completion
criteria are tracked in [TODO.md](TODO.md).

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
