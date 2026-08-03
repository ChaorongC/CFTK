Differential Analysis
=====================

CFTK separates feature-level differential analysis from region-level DMR
analysis.

Feature-Level Analysis
----------------------

Run PCA, differential testing, and visualization for configured modalities:

.. code-block:: bash

   cftk --config cftk_init.json diff

Run one modality:

.. code-block:: bash

   cftk --config cftk_init.json diff --modality cpg

Default matrix locations are derived from ``output_dir`` and the modality name.
For example, ``cpg`` uses:

.. code-block:: text

   <output_dir>/results/1_process/5_merged_matrix/cpg_matrix.tsv

DMR Analysis
------------

Run DMR preparation, ``metilene``, annotation, and plotting:

.. code-block:: bash

   cftk --config cftk_init.json dmr

DMR sample subsets can be configured in ``analysis.dmr.samples``. BedGraph
files are resolved from:

.. code-block:: text

   <output_dir>/results/1_process/4_methylation/

Regenerate Plots
----------------

.. code-block:: bash

   cftk --config cftk_init.json vis --mode diff dmr
