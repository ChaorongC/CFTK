Module Overview
===============

CFTK currently uses script-style modules under ``src``. Public APIs are still
evolving, so this page describes module responsibilities rather than promising
stable Python imports.

Command And Configuration
-------------------------

- ``src/cftk.py``: command-line parser and dispatcher.
- ``src/doctor.py``: non-repairing project, reference, input, and tool
  readiness diagnostics.
- ``src/init.py``: config validation, sample helpers, output path helpers, and
  config summary.
- ``src/process.py``: raw processing steps and CpG matrix merge.
- ``src/run_workflow.py``: fail-fast beginner stage orchestration, artifact
  contracts, resume/quarantine decisions, and per-attempt provenance.
- ``src/util.py``: shared helpers and exact external-command recording.

Analysis Modules
----------------

- ``src/analysis/qc.py``: QC calculations.
- ``src/analysis/power_analysis.py``: statistical power analysis.
- ``src/analysis/model_power.py`` and ``src/analysis/model_power_reference.py``:
  model-development simulation and external aggregate-reference loading.
- ``src/analysis/model_power_discovery.py`` and
  ``src/analysis/model_power_operating_characteristics.py``: cross-validated
  discovery power and null-calibrated operating characteristics.
- ``src/analysis/differential.py`` and ``src/analysis/pca_analysis.py``:
  feature-level differential analysis and PCA.
- ``src/analysis/dmr.py``: DMR preparation, metilene execution, and annotation
  handoff.
- ``src/analysis/occupancy.py`` and ``src/analysis/wps.py``: fragmentomics
  occupancy and WPS features.
- ``src/analysis/delfi.py``, ``src/analysis/end_motif.py``, and
  ``src/analysis/cleavage.py``: finaletoolkit-backed fragmentomics features.
- ``src/analysis/mesa.py`` and ``src/analysis/modality_performance.py``:
  modality performance and MESA modeling.

Visualization And Reporting
---------------------------

- ``src/visualization/``: plotting modules for QC, differential analysis, DMR,
  fragmentomics, MESA, and power analysis.
- ``src/report/report_generator.py``: self-contained HTML report generation.
