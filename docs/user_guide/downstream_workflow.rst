Downstream Workflow
===================

``cftk run`` finishes core processing and QC. Use ``cftk plan`` before any
downstream analysis so CFTK can inspect required matrices, BAMs, reference
components, optional Python packages, external executables, CPU allocation,
and the exact output contract without processing data.

.. code-block:: bash

   cftk plan
   cftk analyze --dry-run
   cftk analyze

``auto`` is intentionally bounded. For a one-group project it selects
occupancy, WPS, and an HTML report. For an explicitly role-defined two-group
project it adds differential analysis. DMR calling and MESA modeling are
available through explicit presets because they introduce R/external-tool and
modeling dependencies that should be inspected before a large run.

When configured differential or MESA modalities include occupancy or WPS and
their matrices do not yet exist, CFTK adds the corresponding feature stage and
runs it before the dependent analysis. The resolved order and dependency edges
are recorded in the plan; other configured modalities must already provide a
matrix or preflight remains blocked.

Targeted-panel scope
--------------------

The default Twist Human Methylome assay is a capture panel. For the WPS,
occupancy, and DELFI stages, ``auto`` therefore creates a deterministic scope
under ``results/4_fragmentomics/_scope/``: reads are filtered to alignments
overlapping the profile ``target_bed``, WPS/occupancy regions are clipped to
panel overlap, and DELFI bins are clipped to panel overlap. The plan, doctor
report, run manifest, stage-local ``fragmentomics_scope.json``, and canonical
``scope.json`` repeat this limitation so panel-level features are not mistaken
for genome-wide measurements. The HTML summary and machine-readable evidence
also show the target identity and derived interval counts after execution. Use
``--fragmentomics-scope genome`` only for validated whole-genome inputs, or
``--fragmentomics-scope panel`` for a custom targeted profile.

For the DMR preset, preflight verifies that ``Rscript`` can load the packages
used by CFTK's bundled annotation script: ``annotatr``, ``GenomicRanges``,
``org.Hs.eg.db``, and ``TxDb.Hsapiens.UCSC.hg38.knownGene``.

.. code-block:: bash

   cftk plan --preset comparative
   cftk analyze --preset dmr
   cftk analyze --preset mesa
   cftk analyze --preset fragmentomics
   cftk analyze --preset all

The ``comparative`` and ``all`` presets require exactly one control group and
one case group from the schema-v2 sample-sheet ``role`` column. CFTK does not
derive labels from group names. One-group occupancy and WPS analyses remain
valid descriptive workflows.

.. figure:: ../_static/cftk_workflow.png
   :alt: CFTK processing, QC, downstream analysis, and report workflow
   :width: 100%

   Downstream analysis starts only after the processed BAMs, CpG bedGraphs,
   matrices, and QC outputs needed by the selected stages are present.

Expected Outputs
----------------

Every ``cftk analyze`` attempt writes an immutable record under
``results/provenance/analysis-runs/<run-id>/``. It contains ``run.json``,
``doctor-before.json``, ``analysis-plan.json``, exact external command records,
``expected-outputs.tsv``, ``figures.tsv``, ``resource-plan.json``, an HTML
summary, and an ``evidence/`` directory with artifact, stage, command, and
resource summaries.

Selected stages add their existing result files without changing their output
locations:

- differential: ``results/3_differential/<modality>/differential_result.tsv``,
  PCA tables, and PCA/violin/heatmap PNG and PDF files;
- DMR: ``results/3_differential/dmr/metilene_input.bedGraph``, raw and
  annotated DMR BED files, and a volcano plot;
- fragmentomics: per-sample occupancy, WPS, DELFI, end-motif, or cleavage
  outputs under ``results/4_fragmentomics/`` with corresponding figures where
  the selected workflow provides them;
- MESA: ``modality_performance.tsv``, ``MESA_model.pkl``,
  ``loocv_predictions.tsv``, and ROC/heatmap/Spearman figures under
  ``results/5_mesa/``;
- report: ``results/report/report.html``.

Optional Per-Sample Jobs
------------------------

The default execution model is local, bounded in-process parallelism through
``--parallel``. For an expensive fragmentomics stage, generate a per-sample
task plan instead:

.. code-block:: bash

   cftk plan --stage delfi end_motif --execution per-sample

Each generated task runs one sample with ``--parallel 1 --no-finalize`` and
therefore has its own terminal status, command-ledger entries, and retry
boundary. After every task for a stage succeeds, its generated finalizer checks
the full cohort, creates any cohort matrix/figures, and adopts the result into
the standard immutable ``cftk analyze`` manifest. Do not run a finalizer after
a failed or incomplete sample task.

The plan is scheduler-neutral: its executable scripts can be submitted by
Slurm, PBS, SGE, or an institutional workflow engine. Clusters using Slurm may
request an array helper:

.. code-block:: bash

   cftk plan --stage delfi --execution per-sample --slurm
   CFTK_SBATCH_ARGS='--cpus-per-task=8 --mem=48G --time=12:00:00' \\
     results/provenance/job-plans/fragmentomics-*/slurm/submit.sh

The helper is written but never submitted automatically. Set its resource
arguments from the chosen stage and local policy. It submits one array element
per sample and a success-gated finalizer job. This does not change FinaleToolkit
parameters or CFTK output schemas; it only moves independent sample work into
separate scheduler jobs.

Resume And Method Boundaries
----------------------------

Successful stages resume when the configuration, lock, and required artifacts
agree. A complete stage recorded by any compatible prior analysis selection is
validated and reused automatically when a later full suite includes it.
Existing untracked outputs with no trusted stage manifest still require
``--adopt-existing``; partial retry outputs are quarantined under
``results/provenance/quarantine/``.

This orchestration layer does not change the scientific implementation of its
stages. The current differential command uses feature-level Mann-Whitney tests
with BH correction. The current DMR stage uses metilene followed by the bundled
R annotation route; it is not the covariate-adjusted beta-regression/DMRcate
workflow. MESA uses the installed ``mesa-cfdna`` package with explicit labels.
Choose or validate a different scientific method explicitly rather than
assuming that a workflow preset changes it.
