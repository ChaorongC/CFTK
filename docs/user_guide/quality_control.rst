cfDNA Quality Control
=====================

The ``qc`` command runs quality-control for three cfDNA features.

1. Methylation distribution
   Uses the merged CpG matrix from raw processing to visualize the methylation distribution in a β-     value density score. Validate the methylation quality in cfDNA samples.

2. Fragment length distribution
   Uses BAM files resolved from mark-duplicate outputs to generate a cfDNA fragent length distribution line plot. Validate the cfDNA fragment length peak range.

3. Dinucleotide frequency
   Requires ``reference_data.genome_fa`` and configured fragment settings.

.. code-block:: bash

   cftk --config cftk_init.json qc -s 1 2 3



Expected Output Location
------------------------

QC outputs are written under:

.. code-block:: text

   <output_dir>/results/2_qc/
