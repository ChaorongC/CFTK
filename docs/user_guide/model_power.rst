Model-Development Power
=======================

CFTK estimates the probability that a fixed biomarker-discovery pipeline
reaches a target out-of-fold cross-validated AUC for a proposed total study
sample size. Filtering, imputation, feature ranking, feature selection,
scaling, and model fitting are repeated inside every training fold.

The ``power`` output is the fraction of simulated studies whose CV AUC reaches
``target_auc``. It evaluates internal model-development adequacy and does not
estimate external generalizability.

This workflow is not the ``cftk power`` command. That legacy CLI command
performs CpG-level analytical power calculations from a configured pickled
standard-deviation table. Model-development power is exposed through the
Python API and Streamlit calculator and uses the manifest-backed NumPy arrays.

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

End-To-End API Smoke Test
-------------------------

The following small calculation verifies reference loading, template
preparation, cross-validation, and matched null calibration. Run it from a
repository checkout after ``python -m pip install .``:

.. code-block:: python

   from analysis.model_power import (
       load_default_model_power_reference,
       prepare_template_ensemble,
   )
   from analysis.model_power_operating_characteristics import (
       run_power_sample_size_grid,
   )

   reference = load_default_model_power_reference(
       reference_dir="data",
       depths=[10],
       sd_stats=["mean"],
       include_index=False,
   )

   templates = prepare_template_ensemble(
       reference.cpg_std_summary,
       reference.cpg_mean,
       n_templates=1,
       template_kwargs={
           "depth": [10],
           "n_features": 20,
           "n_signal_cpgs": 4,
           "meth_diff": 0.06,
           "effect_sd": 0.01,
       },
       random_state=42,
   )

   result = run_power_sample_size_grid(
       templates,
       sample_sizes=[20],
       simulations_per_template=1,
       null_simulations_per_template=19,
       power_kwargs={
           "models": ("logreg",),
           "cv_folds": 2,
           "top_k": 3,
           "target_auc": 0.70,
       },
       alpha=0.05,
       ci_method="none",
       n_jobs=1,
       random_state=43,
   )

   columns = [
       "sample_size",
       "mean_depth",
       "mean_cv_auc",
       "detection_power",
       "target_attainment_probability",
       "probability_of_success",
   ]
   print(result["power_curve"][columns])

The result also contains ``replicate_results``, ``null_replicate_results``,
``template_summary``, and ``run_metadata``. This one-template, one-simulation
example is only a wiring smoke test. Increase template and simulation counts
for study-design analysis; use substantially more null simulations when precise
tail probabilities are required.

Use ``ci_method='none'`` for fast web calculations. Pooled Wilson and
hierarchical bootstrap intervals are available when uncertainty estimates are
required.
