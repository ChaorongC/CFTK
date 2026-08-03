Model-Development Power
=======================

CFTK estimates the probability that a fixed biomarker-discovery pipeline
reaches a target out-of-fold cross-validated AUC for a proposed total study
sample size. Filtering, imputation, feature ranking, feature selection,
scaling, and model fitting are repeated inside every training fold.

The ``power`` output is the fraction of simulated studies whose CV AUC reaches
``target_auc``. It evaluates internal model-development adequacy and does not
estimate external generalizability.

For an interactive interface to this workflow, see
:doc:`model_power_calculator`.

The main functions are:

- ``analysis.model_power.load_default_model_power_reference``
- ``analysis.model_power.prepare_template_ensemble``
- ``analysis.model_power_discovery.run_power_sample_size_grid``
- ``analysis.model_power_operating_characteristics.run_power_sample_size_grid``
- ``visualization.plot_model_power.plot_power_by_sample_size``

The Python distribution contains these APIs but not the aggregate reference
arrays. Load repository-hosted or separately obtained arrays explicitly:

.. code-block:: python

   from analysis.model_power import load_default_model_power_reference

   reference = load_default_model_power_reference(
       reference_dir="/path/to/CFTK/data",
       depths=[10, 30],
       sd_stats=["mean"],
       include_index=False,
   )

``CFTK_MODEL_POWER_DATA`` can provide the same directory for installed callers.
The GitHub Streamlit application passes its repository-local ``data/`` path.

Use ``ci_method='none'`` for fast web calculations. Pooled Wilson and
hierarchical bootstrap intervals are available when uncertainty estimates are
required.
