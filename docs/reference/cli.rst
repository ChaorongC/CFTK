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
   Create a missing schema-v2 project interactively, or validate an existing
   schema-v2/legacy config and prepare bwa-meth, samtools, and Picard reference
   companions.

   .. code-block:: bash

      cftk init

   For managed batch setup, pass ``--non-interactive`` and
   ``--sample-sheet PATH``. The default profile installs under
   ``CFTK_REFERENCE_ROOT`` or ``~/.cache/cftk/references``. Local mode also
   requires ``--reference-mode local --reference-root PATH``. ``--profile`` and
   ``--profile-version`` select a non-default or ambiguous version. The expert
   ``CFTK_REFERENCE_REGISTRY`` override remains available. Pass
   ``--skip-reference-prep`` for validation only.

``doctor``
   Check whether selected core processing and optional downstream-analysis
   stages can start. The default checks steps 1 through 4. Diagnostics continue
   after failures and do not download, index, repair, or process data.

   .. code-block:: bash

      cftk --config cftk_init.json doctor
      cftk --config cftk_init.json doctor --step 4 --json
      cftk --config cftk_init.json doctor --analysis-preset comparative

   ``--step {1,2,3,4}`` accepts one or more space-separated values.
   ``--target-bed PATH`` validates an expert Picard covered-target override;
   ``--skip-picard-metrics`` removes Picard metrics from the required step-3
   tool checks. ``--parallel N`` validates a parallel-sample override against
   the configured total CPU budget and detected scheduler allocation.
   ``--json`` writes only machine-readable JSON to stdout.
   ``--analysis-preset`` or ``--analysis-stage`` adds read-only checks for the
   downstream stage dependencies, inputs, references, roles, and CPU plan.
   ``--fragmentomics-scope`` applies the same assay-aware scope resolution to
   selected WPS, occupancy, and DELFI checks.

   Human checks use ``PASS``, ``WARN``, and ``FAIL``. No required failures
   returns exit status 0; one or more failures returns 1; invalid arguments
   return 2.

``run``
   Run the beginner-safe schema-v2 core workflow. The default order is process
   steps 1-4, then QC steps 2, 0, and 1. FASTQ projects run every stage; BAM
   projects record process steps 1-2 as skipped. Mixed inputs, legacy configs,
   missing locks, alternative processing tools, and failed doctor checks stop
   before computation.

   .. code-block:: bash

      cftk run
      cftk run --dry-run
      cftk run --parallel 2 --qc-dinucleotide

   ``--adopt-existing`` is required to validate complete outputs that predate a
   run manifest. Partial untrusted outputs are preserved in a timestamped
   quarantine before retry. Automatic resume requires a matching config, lock,
   options hash, trusted prior stage state, and currently valid artifacts.
   ``--target-bed PATH`` is an expert one-run Picard target override.
   ``--parallel N`` sets concurrent samples. CFTK divides the configured total
   ``process.cores`` budget across those sample commands and records the result
   in ``resource-plan.json``.

   Every attempt writes ``run.json``, event and command JSONL files, doctor and
   tool-version JSON, output/figure TSVs, and ``run-summary.html`` under
   ``results/provenance/runs/<run-id>/``. It also generates an ``evidence/``
   directory containing stage/artifact/command tables and sanitized figure
   previews. If analysis completes but evidence reporting fails, the manifest
   status is ``complete_with_reporting_error`` and the command returns a
   distinct nonzero status; a subsequent run can rebuild evidence without
   rerunning valid stages. See
   :doc:`../user_guide/beginner_run` for the full contract.

``plan``
   Resolve a role-aware downstream preset and record its dependencies,
   resources, expected outputs, and read-only doctor result without launching
   a stage.

   .. code-block:: bash

      cftk plan
      cftk plan --preset all
      cftk plan --stage diff report --json

   ``auto`` plans occupancy, WPS, and reporting for one group; it additionally
   plans differential analysis for an explicit two-group control/case project.
   ``comparative`` and ``all`` require roles, not group-name inference. Plans
   are recorded under ``results/provenance/analysis-plans/``. Use
   ``--fragmentomics-scope panel|genome`` only when overriding the assay-aware
   default. If configured
   differential or MESA modalities require occupancy or WPS matrices, their
   producer stages are added ahead of the dependent stage.

``analyze``
   Run downstream stages with fail-fast preflight, artifact contracts,
   provenance, evidence, and resume behavior. It requires a schema-v2 project
   and its matching lock file.

   .. code-block:: bash

      cftk analyze --dry-run
      cftk analyze --preset fragmentomics
      cftk analyze --preset comparative
      cftk analyze --preset all

   Explicit presets are ``descriptive``, ``differential``, ``dmr``,
   ``fragmentomics``, ``mesa``, ``comparative``, ``all``, and ``report``.
   Use ``--stage`` for a precise stage or alias such as ``diff``, ``wps``, or
   ``report``. ``--fragmentomics-scope`` controls targeted WPS, occupancy, and
   DELFI inputs; it defaults to panel scope for Twist. ``--adopt-existing`` validates complete outputs produced before
   an analysis manifest; untrusted partial outputs are quarantined before a
   retry. See :doc:`../user_guide/downstream_workflow`.

``process``
   Run raw processing steps 1 through 4. Step 3 uses the schema-v2 profile's
   covered-target BED for Picard metrics, with a bundled source-checkout
   fallback for legacy configs.

   .. code-block:: bash

      cftk --config cftk_init.json process -s 1 2 3 4

   Use ``--target-bed PATH`` to override the covered targets, or
   ``--skip-picard-metrics`` for a workflow that does not need them.

``qc``
   Run QC steps 0 through 3. Step 0 assembles QC tables, step 1 plots
   methylation distributions, step 2 measures fragment lengths, and step 3
   computes dinucleotide frequencies.

   .. code-block:: bash

      cftk --config cftk_init.json qc -s 1 2 3

``power``
   Run the legacy CpG-level analytical power workflow. It uses the pickled
   ``reference_data.cpg_std`` table from ``cftk_init.json`` to evaluate
   methylation effect-size and sample-size scenarios.

   .. code-block:: bash

      cftk --config cftk_init.json power -s 100 -e 0.1

   This command is separate from the model-development power calculator. The
   latter simulates cross-validated feature selection and classification,
   performs matched null calibration, and is available through the Python API
   and Streamlit app described in :doc:`../user_guide/model_power`.

``diff``
   Run PCA, differential testing, and differential visualizations.

   .. code-block:: bash

      cftk --config cftk_init.json diff --modality cpg

``dmr``
   Run DMR analysis.

   .. code-block:: bash

      cftk --config cftk_init.json dmr

``frag``
   Run fragmentomics workflows. Occupancy and WPS are valid for a one-group
   descriptive project; comparison figures are produced only when two groups
   are available. The default Twist profile scopes WPS, occupancy, and DELFI to
   panel-overlapping reads and regions. Use ``--fragmentomics-scope genome``
   only for validated whole-genome inputs.

   .. code-block:: bash

      cftk --config cftk_init.json frag --wps
      cftk --config cftk_init.json frag --delfi --fragmentomics-scope panel

``mesa``
   Run MESA modality performance, model construction, and LOOCV. It requires
   explicit control/case roles; CFTK does not infer labels from group names.

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
   Run the expert compatibility pipeline. It may continue after a failed step
   and does not provide the validated run-state/resume contract of ``cftk run``.

   .. code-block:: bash

      cftk --config cftk_init.json run-all --parallel 4
