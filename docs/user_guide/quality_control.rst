cfDNA Quality Control
=====================

The ``qc`` command has one table-assembly step and three cfDNA feature checks.

0. QC metrics and scores
   Parses process reports into ``qc_summary.tsv`` and ``qc_scores.tsv``. Run
   fragment-length step 2 first when the summary should include median fragment
   length.

1. Methylation distribution
   Uses the merged CpG matrix to plot cohort methylation beta-value densities.
   Interpret this with coverage, M-bias, conversion, and assay information; its
   shape alone is not a sample-quality verdict.

2. Fragment length distribution
   Uses duplicate-marked BAMs and deepTools ``bamPEFragmentSize`` to write raw
   lengths, per-sample histograms, and a combined cohort plot.

3. Dinucleotide frequency
   Requires the reference FASTA and configured fragment settings. This
   expensive stage is opt-in for ``cftk run``.

.. code-block:: bash

   cftk --config cftk_init.json qc -s 2 0 1

The beginner workflow uses that order by default. Add step 3 explicitly when
needed:

.. code-block:: bash

   cftk --config cftk_init.json qc -s 3



Expected Output Location
------------------------

QC outputs are written under:

.. code-block:: text

   <output_dir>/results/2_qc/
