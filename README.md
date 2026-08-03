# CFTK

CFTK is a cfDNA multimodal epigenetic analysis toolkit for processing
cfMethyl-Seq style data and running downstream methylation, fragmentomics,
visualization, modeling, and report workflows. 
For detailed guidance and tutorial, please refer to the [CFTK website](https://chaorongc.github.io/CFTK/index.html)

You can also use the [CFTK model power calculator](https://cftk-model-power.streamlit.app/) before you start to process your cfDNA cohort.


The package is under active development. The current command-line entry point is
implemented in `src/cftk.py` and is driven by a project configuration file named
`cftk_init.json`.

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

Validate the example configuration:

```bash
cftk --config cftk_init.json init
```

Inspect available commands:

```bash
cftk --help
```

Run a raw processing step after editing `cftk_init.json` for your samples,
reference files, tools, and output directory:

```bash
cftk --config cftk_init.json process -s 1 2 3 4
```

Some workflows require external bioinformatics tools and reference files that
are not installed by Python packaging alone. See the documentation for details.

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
