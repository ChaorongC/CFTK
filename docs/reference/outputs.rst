Output Layout
=============

CFTK writes results below ``<output_dir>/results``.

.. list-table::
   :header-rows: 1

   * - Path
     - Purpose
   * - ``0_power/``
     - Statistical power analysis outputs and plots.
   * - ``provenance/commands.jsonl``
     - Append-only start and finish records for exact external commands,
       including return codes and working directories.
   * - ``provenance/runs/<run-id>/``
     - Immutable beginner-run attempt directory containing stage state,
       config/lock/options identity, doctor and tool reports, an immediate
       command-ledger mirror, expected output/figure tables, events, and an
       HTML summary.
   * - ``provenance/runs/<run-id>/resource-plan.json``
     - Resolved total CPU budgets, parallel samples, per-sample threads,
       estimated peak threads, and detected scheduler allocation by stage.
   * - ``provenance/runs/<run-id>/evidence/``
     - Automatically generated stage/artifact/command TSVs and evidence
       figures for every terminal ``cftk run`` attempt. Tables may contain
       private absolute paths and verbatim commands; do not commit them. The
       source-checkout helper remains available for historical manifests.
   * - ``provenance/runs/<run-id>/run-summary.html``
     - Core processing/QC summary. When ``cftk run --downstream PRESET`` is
       used, it also links the downstream manifest and HTML summary and records
       the downstream preset/status without merging the two provenance files.
   * - ``provenance/analysis-plans/<plan-id>/``
     - Read-only downstream dependency, input, role, resource, and output
       plan written by ``cftk plan``.
   * - ``provenance/analysis-runs/<run-id>/``
     - Immutable downstream-analysis attempt containing the selected preset,
       doctor preflight, artifact/figure contracts, command mirror, resource
       plan, evidence bundle, and HTML summary.
   * - ``provenance/job-plans/<plan-id>/``
     - Advanced per-sample task scripts, finalizer scripts, and ``job-plan.json``
       written by ``cftk plan --execution per-sample``. They are never
       submitted automatically.
   * - ``provenance/job-finalizers/<plan-id>/``
     - Completion markers written only after the generated cohort finalizer
       succeeds. ``cftk status`` checks these markers together with current
       sample artifacts; it does not replace scheduler status commands.
   * - ``provenance/quarantine/<run-id>/``
     - Preserved partial or damaged stage artifacts moved before a safe retry.
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
     - QC result tables and figures. The beginner guide includes rendered
       examples for the fragment-length and methylation-distribution outputs.
   * - ``3_differential/``
     - Differential result tables, PCA intermediates/figures, heatmaps,
       violin plots, and DMR outputs. Managed differential manifests also
       record the selected modalities and SHA-256 input-matrix signatures.
   * - ``4_fragmentomics/occupancy/``
     - Occupancy features and ``fragmentomics_scope.json`` when scoped.
   * - ``4_fragmentomics/wps/``
     - WPS features and ``fragmentomics_scope.json`` when scoped.
   * - ``4_fragmentomics/delfi/``
     - DELFI-style features and ``fragmentomics_scope.json`` when scoped.
   * - ``4_fragmentomics/end_motif/``
     - End motif features.
   * - ``4_fragmentomics/cleavage/``
     - Cleavage features.
   * - ``4_fragmentomics/_scope/<scope-id>/``
     - Deterministic targeted-fragmentomics intermediates: clipped panel
       regions/bins, panel-read BAMs and indexes, and ``scope.json``. These
       files document the panel limitation and are required for resume.
   * - ``5_mesa/``
     - MESA performance, model, and LOOCV outputs.
   * - ``report/``
     - Self-contained HTML report.

Completion Checks
-----------------

Launching a command is not enough to mark a workflow complete. Check command
exit status, the command provenance ledger, logs, and expected files under the
output directory before using downstream steps.

``cftk run`` automates the stage-level subset of these checks. Its
``expected-outputs.tsv`` and ``figures.tsv`` enumerate the exact required paths,
and ``run.json`` distinguishes ``planned``, ``pending``, ``running``,
``complete``, ``complete_with_reporting_error``, ``failed``, ``interrupted``,
``skipped``, ``resumed``, and ``adopted`` states. Artifact validation checks
required file existence and nonempty content; it does not currently calculate
checksums for large outputs. Differential stage reuse is stricter: every
selected input matrix must match the SHA-256 signature stored in the trusted
analysis manifest. The final report links each differential TSV and shows
modality context, PCA/violin/heatmap figures, and a lowest-q navigation table
without applying a significance threshold.
