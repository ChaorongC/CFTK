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
You need a local profile under one reference root; managed download is not yet
enabled.

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
- ``run-all``: run the configured end-to-end workflow.

4. Run Rawdata Processing
-------------------------

After initialization:

.. code-block:: bash

   cftk --config cftk_init.json process -s 1 2 3 4

The process command creates standard subdirectories under
``<output_dir>/results/1_process`` and merges per-sample CpG calls into
``cpg_matrix.tsv`` after successful methylation calling.

5. Run Downstream Analysis
---------------------------

Examples:

.. code-block:: bash

   cftk --config cftk_init.json qc -s 1 2 3
   cftk --config cftk_init.json diff
   cftk --config cftk_init.json frag --wps
   cftk --config cftk_init.json mesa --performance --mesa-model --loocv
   cftk --config cftk_init.json report


6. End-To-End Runs
------------------

The ``run-all`` command runs the configured pipeline end to end:

.. code-block:: bash

   cftk --config cftk_init.json run-all

Because ``run-all`` continues after some failures, review logs and expected
artifacts before treating a run as complete.
