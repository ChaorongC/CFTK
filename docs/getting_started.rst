Getting Started
===============

CFTK workflow is organized by a project configuration file and worked as a command-line
dispatcher. Edit ``cftk_init.json`` for your samples, then run the commands that
match the analysis stage you need.

1. Json Config Example
==========================

Below is the standard CFTK configuration file used to set the workflow:

.. literalinclude:: _static/cftk_init.json
   :language: json
   :linenos:
   :caption: cftk_init.json
   :class: long-json-block

.. raw:: html

   <style>
   .long-json-block .highlight {
       max-height: 400px !important;
       overflow-y: auto !important;
   }
   </style>


2. Validate The Config
----------------------

.. code-block:: bash

   cftk --config cftk_init.json init

This checks required top-level sections, sample definitions, comparison group
names, and required reference fields.

3. The Help Commands
--------------------

.. code-block:: bash

   cftk --help

The major commands are:

- ``init``: validate and summarize ``cftk_init.json``.
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

After set up all the related paths in ``cftk_init.json``:

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
