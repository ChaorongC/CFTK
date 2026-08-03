Fragmentomics
=============

The ``frag`` command runs fragmentomics workflows. If no sub-workflow flag is
provided, CFTK attempts all configured fragmentomics analyses.

.. code-block:: bash

   cftk --config cftk_init.json frag

Sub-Workflows
-------------

``--occupancy``
   Run DANPOS-style nucleosome occupancy analysis.

``--wps``
   Compute window protection score features.

``--delfi``
   Run DELFI-style fragment ratio features through ``finaletoolkit``.

``--end-motif``
   Run k-mer end motif analysis through ``finaletoolkit``.

``--cleavage``
   Run CTCF cleavage analysis through ``finaletoolkit``.

Examples
--------

Run only WPS:

.. code-block:: bash

   cftk --config cftk_init.json frag --wps

Run occupancy and DELFI:

.. code-block:: bash

   cftk --config cftk_init.json frag --occupancy --delfi

Reference Inputs
----------------

Fragmentomics workflows use different reference files:

- ``chrom_sizes`` for genomic intervals and bigWig/binned workflows.
- ``genome_2bit`` for DELFI and some finaletoolkit commands.
- ``tss_pas_bed`` for WPS and occupancy regions.
- ``ctcf_bed`` for cleavage.
- ``blacklist``, ``gap``, and ``bins`` for DELFI-style features.

Outputs are written under:

.. code-block:: text

   <output_dir>/results/4_fragmentomics/
