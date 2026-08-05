Getting Started
===============

CFTK uses a compact project JSON, a TSV sample sheet, and one versioned
reference profile. Existing legacy nested JSON remains supported.

1. Initialize A Project
=======================

From the project directory, run:

.. code-block:: bash

   cftk init

The guided setup uses current-directory defaults. It may create ``samples.tsv``
from unambiguous FASTQ/BAM inputs and stop so you can assign explicit
``control`` and ``case`` roles. Rerun the same command after editing the sheet.
CFTK selects managed mode automatically and installs the pinned
``twist_human_methylome_hg38`` profile version ``1.0.0`` under
``CFTK_REFERENCE_ROOT`` or ``~/.cache/cftk/references``. Individual reference
paths are not added to the project JSON.

Before starting the workflow, run:

.. code-block:: bash

   cftk doctor

Resolve every ``FAIL`` before processing. ``WARN`` identifies optional tooling
or provenance that deserves review but does not prevent execution. The check
does not repair the project; rerun ``cftk init`` when its remedy requests
reference preparation or lock regeneration.

2. Inspect The Compact Files
----------------------------

See :doc:`user_guide/configuration` for the complete contract. Minimal examples
are stored under ``examples/`` in the repository:

- ``cftk_init.schema-v2.json``
- ``samples.tsv``
- ``reference-profile-manifest.json``

Initialization validates file hashes, target BED coordinates, chromosome
sizes, and the prepared FASTA index. It writes ``cftk.lock.json`` without
machine-specific reference paths.

3. The Help Commands
--------------------

.. code-block:: bash

   cftk --help

The major commands are:

- ``init``: validate the config and prepare the reference genome.
- ``doctor``: check selected processing readiness without changing the project.
- ``run``: run the fail-fast schema-v2 core processing and QC workflow.
- ``process``: run raw processing steps 1 through 4.
- ``qc``: run methylation, fragment length, or dinucleotide QC.
- ``power``: run statistical power analysis.
- ``diff``: run PCA, differential testing, and summary plots.
- ``dmr``: run DMR preparation, metilene, annotation, and plotting.
- ``frag``: run occupancy, WPS, DELFI, end motif, and cleavage workflows.
- ``mesa``: run modality performance and multimodal MESA modeling.
- ``merge``: build feature matrices from user-provided files.
- ``vis``: regenerate plots from existing results.
- ``report``: generate a self-contained HTML report.
- ``run-all``: run the expert compatibility workflow, which may continue after
  failures.

4. Run The Beginner Workflow
----------------------------

After initialization:

.. code-block:: bash

   cftk run

This runs processing steps 1 through 4, fragment-length QC, QC table assembly,
and methylation-distribution QC. It runs doctor first, stops before downstream
stages after a required failure, and validates all required files. BAM-only
projects skip trimming and alignment. Mixed FASTQ/BAM projects are rejected.

Review the human summary at
``results/provenance/runs/<run-id>/run-summary.html``. See
:doc:`user_guide/beginner_run` for dry-run, resume, adoption, quarantine,
Slurm, and exact artifact behavior.

Exact external commands are appended to
``<output_dir>/results/provenance/commands.jsonl`` before execution, followed
by completion records containing exit status. Keep this ledger with the
project configuration, lock file, and outputs when archiving an analysis.

5. Run Expert Or Downstream Commands
------------------------------------

Examples:

Install ``.[analysis]`` and/or ``.[fragmentomics]`` before invoking the
corresponding optional commands shown below.

.. code-block:: bash

   cftk --config cftk_init.json process -s 1 2 3 4
   cftk --config cftk_init.json qc -s 1 2 3
   cftk --config cftk_init.json diff
   cftk --config cftk_init.json frag --wps
   cftk --config cftk_init.json mesa --performance --mesa-model --loocv
   cftk --config cftk_init.json report


6. Expert Compatibility Runs
----------------------------

The legacy ``run-all`` command spans additional configured analyses:

.. code-block:: bash

   cftk --config cftk_init.json run-all

``run-all`` catches some failures and continues, and it does not use the
beginner run-state contract. Keep it for advanced compatibility workflows; use
``cftk run`` when fail-fast completion and validated resume are required.
