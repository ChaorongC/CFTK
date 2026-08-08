Installation
============

CFTK can be installed from a source checkout. External bioinformatics tools
and large reference datasets remain separate from the Python distribution.

A source installation has two required steps:

1. Create and activate the environment. This installs the external
   bioinformatics tools used by the processing workflow.
2. Install CFTK inside that environment. This installs the Python package and
   the ``cftk`` command.

Run both steps in order. Creating the environment alone does not install CFTK.

Clone The Repository
--------------------

.. code-block:: bash

   git clone https://github.com/ChaorongC/CFTK.git
   cd CFTK

Create An Environment
---------------------

The repository includes an environment that pins the core processing
executables and Java runtime used for validation. Install Mamba or Conda, then
create it from the repository root:

.. code-block:: bash

   mamba env create -f environment.yml
   mamba activate cftk

The environment pins Python, Java, Trim Galore, FastQC, bwa-meth, BWA,
Sambamba, samtools, Picard, MethylDackel, BEDTools, deepTools, and MultiQC. It is a
portable environment specification, not a byte-for-byte solver lock. Record an
explicit platform lock or export for a production analysis.

Install CFTK
------------

Install the beginner processing package from the checkout:

.. code-block:: bash

   python -m pip install .

This default installs the Python libraries used by ``cftk init``, ``doctor``,
``run``, process/QC, reporting, and their standard figures. Install optional
workflow groups only when they are needed:

.. code-block:: bash

   python -m pip install ".[analysis]"       # differential, power, and MESA
   python -m pip install ".[fragmentomics]" # WPS and FinaleToolkit workflows
   python -m pip install ".[analysis,fragmentomics]" # expert run-all inputs

Include the calculator's statistical and Streamlit dependencies when
developing or deploying it:

.. code-block:: bash

   python -m pip install ".[web]"

The installation provides the ``cftk`` console command and the model-power
Python modules. It does not install repository-root reference data.

Planned Distribution
--------------------

The two-step source-checkout workflow above is the supported installation path
today. No intermediate change to make ``environment.yml`` install the local
checkout automatically is planned or required.

The eventual distribution target is a versioned Conda package, preferably
through Bioconda. It will install CFTK and its validated default processing
toolchain with one Conda or Mamba command and no repository checkout. Analysis,
fragmentomics, and web dependencies will remain optional rather than expanding
the default beginner environment. See the distribution milestones in the
repository ``TODO.md``.

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
   mkdir example_study
   cd example_study
   cftk init
   cftk doctor
   cftk run --dry-run
   cftk run

The guided initializer creates compact project metadata and installs the pinned
managed default profile. Acquisition peaks at approximately 5.7 GB before
temporary download artifacts are removed; bwa-meth indexing requires additional
space. CFTK tries ``bwameth index`` first and ``bwameth.py index`` if that fails,
then prepares ``.fai`` and Picard ``.dict`` companions. Use
``--skip-reference-prep`` only when those files are managed outside CFTK.

``cftk doctor`` is read-only apart from an ephemeral output-location write
probe. It verifies the selected process tools, full profile checksums, reference
companions, project lock, FASTQ/BAM inputs, BAM/reference sequence dictionaries,
and output capacity. It never downloads references, creates indexes, repairs
files, or changes source data. Full checksum verification reads several
gigabytes and can take a few minutes on shared storage.

For a downstream BAM-only project, check only methylation-calling readiness:

.. code-block:: bash

   cftk doctor --step 4

Machine-readable output is available for schedulers:

.. code-block:: bash

   cftk doctor --json > doctor.json

The command exits 0 when there are no required failures, including reports
that contain only optional warnings. It exits 1 for readiness failures;
argument errors retain argparse exit status 2.

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
- Picard (the ``picard`` executable)
- ``MethylDackel``
- ``bedtools``
- ``multiqc``
- ``bamPEFragmentSize`` from deepTools
- ``DANPOS``
- UCSC tools such as ``wigToBigWig`` and ``bigWigAverageOverBed``
- ``metilene`` for DMR analysis
- R packages used by DMR annotation, including ``annotatr``, ``GenomicRanges``,
  ``org.Hs.eg.db``, and the hg38 ``TxDb`` annotation package

The pinned environment installs the core step 1-4 tools. Install advanced
workflow tools separately in the compute environment where those workflows
will run.

Build The Documentation
-----------------------

.. code-block:: bash

   python -m pip install -r docs/requirements.txt
   python -m sphinx -W -b html docs docs/_build/html

The local HTML entry point is ``docs/_build/html/index.html``.
