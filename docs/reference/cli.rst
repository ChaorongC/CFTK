Command Reference
=================

Global Options
--------------

.. code-block:: text

   cftk [--config PATH] <command> ...

``--config PATH``
   Path to ``cftk_init.json``. Defaults to ``./cftk_init.json``.

Commands
--------

``init``
   Validate ``cftk_init.json`` and print a project summary.

   .. code-block:: bash

      cftk --config cftk_init.json init

``process``
   Run raw processing steps 1 through 4.

   .. code-block:: bash

      cftk --config cftk_init.json process -s 1 2 3 4

``qc``
   Run QC steps 1 through 3.

   .. code-block:: bash

      cftk --config cftk_init.json qc -s 1 2 3

``power``
   Run statistical power analysis.

   .. code-block:: bash

      cftk --config cftk_init.json power -s 100 -e 0.1

``diff``
   Run PCA, differential testing, and differential visualizations.

   .. code-block:: bash

      cftk --config cftk_init.json diff --modality cpg

``dmr``
   Run DMR analysis.

   .. code-block:: bash

      cftk --config cftk_init.json dmr

``frag``
   Run fragmentomics workflows.

   .. code-block:: bash

      cftk --config cftk_init.json frag --wps

``mesa``
   Run MESA modality performance, model construction, and LOOCV.

   .. code-block:: bash

      cftk --config cftk_init.json mesa --performance --mesa-model --loocv

``merge``
   Build feature matrices from user-specified files in the config ``merge``
   block.

   .. code-block:: bash

      cftk --config cftk_init.json merge --modality cpg

``vis``
   Regenerate visualizations from existing results.

   .. code-block:: bash

      cftk --config cftk_init.json vis --mode all

``report``
   Generate a self-contained HTML report.

   .. code-block:: bash

      cftk --config cftk_init.json report

``run-all``
   Run the configured end-to-end pipeline.

   .. code-block:: bash

      cftk --config cftk_init.json run-all --parallel 4
