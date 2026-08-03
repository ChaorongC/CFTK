Raw Processing
==============

The ``process`` command runs raw processing steps 1 through 4.

.. code-block:: bash

   cftk --config cftk_init.json process -s 1 2 3 4

Steps
-----

1. Adapter trimming
   Uses ``trim_galore`` or ``fastp`` for FASTQ inputs.

2. Bisulfite alignment
   Uses ``bwameth`` or ``bismark``.

3. Duplicate marking
   Uses ``sambamba``, ``picard``, or ``samblaster``.

4. CpG methylation calling
   Uses ``MethylDackel`` or ``bismark_methylation_extractor``.

Parallel Samples
----------------

Use ``--parallel`` to process multiple samples concurrently per step:

.. code-block:: bash

   cftk --config cftk_init.json process -s 1 2 3 4 --parallel 4

CFTK splits configured cores across parallel samples. For example, if a step
uses 20 total cores and ``--parallel 4`` is set, each sample receives 5 cores.

Merged CpG Matrix
-----------------

After step 4, CFTK can merge per-sample CpG bedGraph files into:

.. code-block:: text

   <output_dir>/results/1_process/5_merged_matrix/cpg_matrix.tsv

The merged matrix is the default input for methylation QC, differential
analysis, MESA modeling, and report generation.

Validation Strategy
-------------------

For a new compute environment, validate steps incrementally:

.. code-block:: bash

   cftk --config cftk_init.json process -s 1
   cftk --config cftk_init.json process -s 2
   cftk --config cftk_init.json process -s 3
   cftk --config cftk_init.json process -s 4

Check logs and expected output files after each step before using ``run-all``.
