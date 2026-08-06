User Guide
==========

The user guide follows the order of a typical CFTK project: configure samples
and references, process raw data, run quality control, compute fragmentomics
features, perform differential analysis, model modalities, and generate reports.

Expected Outputs
----------------

Each page below embeds a visible example and names the files or in-memory
objects that a real run produces. The processing/QC examples marked as
sanitized are observed technical workflow evidence. Downstream examples marked
synthetic are fixed-seed illustrations for learning the output shape only.

.. figure:: ../_static/cftk_workflow.png
   :alt: CFTK workflow from inputs through processing, QC, analysis, and report
   :width: 100%

   High-level orientation: the result roots are ``results/1_process``,
   ``results/2_qc``, ``results/3_differential``,
   ``results/4_fragmentomics``, ``results/5_mesa``, and ``results/report``.

.. toctree::
   :maxdepth: 2

   configuration
   beginner_run
   processing
   quality_control
   differential
   fragmentomics
   MESA_modeling
   model_power
   model_power_calculator
   validation_acceptance
   visualization_report
