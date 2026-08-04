Output Layout
=============

CFTK writes results below ``<output_dir>/results``.

.. list-table::
   :header-rows: 1

   * - Path
     - Purpose
   * - ``0_power/``
     - Statistical power analysis outputs and plots.
   * - ``1_process/1_trimming/``
     - Trimmed FASTQ files and trimming QC reports.
   * - ``1_process/2_alignment/``
     - Aligned BAM files and alignment metrics.
   * - ``1_process/3_markdup/``
     - Duplicate-marked BAM files, indexes, and ``picard_metrics/`` target and
       alignment performance outputs.
   * - ``1_process/4_methylation/``
     - Per-sample CpG methylation calls.
   * - ``1_process/5_merged_matrix/``
     - Merged CpG matrix.
   * - ``2_qc/``
     - QC result tables and figures.
   * - ``3_differential/``
     - Differential, PCA, heatmap, violin, and DMR outputs.
   * - ``4_fragmentomics/occupancy/``
     - Occupancy features.
   * - ``4_fragmentomics/wps/``
     - WPS features.
   * - ``4_fragmentomics/delfi/``
     - DELFI-style features.
   * - ``4_fragmentomics/end_motif/``
     - End motif features.
   * - ``4_fragmentomics/cleavage/``
     - Cleavage features.
   * - ``5_mesa/``
     - MESA performance, model, and LOOCV outputs.
   * - ``report/``
     - Self-contained HTML report.

Completion Checks
-----------------

Launching a command is not enough to mark a workflow complete. Check command
exit status, logs, and expected files under the output directory before using
downstream steps.
