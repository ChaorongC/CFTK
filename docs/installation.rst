Installation
============

CFTK can be installed from a source checkout. External bioinformatics tools
and large reference datasets remain separate from the Python distribution.

Clone The Repository
--------------------

.. code-block:: bash

   git clone https://github.com/ChaorongC/CFTK.git
   cd CFTK

Create An Environment
---------------------

Use a Python environment that matches the target platform for your analysis.
The project metadata currently declares Python 3.9 or newer.

.. code-block:: bash

   python -m venv .venv
   source .venv/bin/activate
   python -m pip install --upgrade pip

Install Python dependencies according to the workflows you plan to run. The
core analysis modules use packages such as ``numpy``, ``pandas``, ``scipy``,
``scikit-learn``, ``matplotlib``, ``seaborn``, ``pysam``, ``pyBigWig``,
``bx-python``, ``statsmodels``, ``xgboost``, and ``finaletoolkit``.

Install CFTK
------------

Install the core package from the checkout:

.. code-block:: bash

   python -m pip install .

Include the Streamlit dependency when developing or deploying the calculator:

.. code-block:: bash

   python -m pip install ".[web]"

The installation provides the ``cftk`` console command and the model-power
Python modules. It does not install repository-root reference data.

Run The Model-Power Calculator
------------------------------

From the repository root, install the web extra and start Streamlit:

.. code-block:: bash

   python -m pip install ".[web]"
   streamlit run apps/model_power_calculator.py

The app resolves ``data/`` relative to the repository root. Keep the app and
reference-data directories in their checkout layout. Streamlit prints the
local browser URL after startup, normally ``http://localhost:8501``.

Run CFTK
--------

.. code-block:: bash

   cftk --help
   cftk --config cftk_init.json init

Direct source execution remains supported:

.. code-block:: bash

   python src/cftk.py --help

Model-Power Reference Data
--------------------------

The calculator algorithms and loaders are installed, but the aggregate CpG
arrays under ``data/`` are intentionally excluded from wheels and source
distributions. The 24 tracked calculator arrays occupy about 520 MB (496 MiB).
Supply a repository checkout explicitly:

.. code-block:: python

   from analysis.model_power import load_default_model_power_reference

   reference = load_default_model_power_reference(
       reference_dir="/path/to/CFTK/data",
       depths=[10, 30],
       sd_stats=["mean"],
       include_index=False,
   )

Alternatively, configure the directory once for a process:

.. code-block:: bash

   export CFTK_MODEL_POWER_DATA=/path/to/CFTK/data

An explicit ``reference_dir`` takes precedence over the environment variable.
The GitHub Streamlit app passes its checked-out ``data/`` directory directly.

External Tools
--------------

Many workflows call command-line tools that Python packaging does not install:

- ``trim_galore``
- ``bwameth``
- ``sambamba``, ``samtools``
- ``MethylDackel``
- ``bedtools``
- ``multiqc``
- ``DANPOS``
- UCSC tools such as ``wigToBigWig`` and ``bigWigAverageOverBed``
- ``metilene`` for DMR analysis
- R packages used by DMR annotation, including ``annotatr`` and hg38 annotation
  packages

Install and validate these tools separately in the compute environment where
the pipeline will run.

Build The Documentation
-----------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python -m sphinx -W -b html docs docs/_build/html

The local HTML entry point is ``docs/_build/html/index.html``.
