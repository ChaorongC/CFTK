Reference Data
==============

CFTK needs both repository-hosted region resources and user-supplied genome
resources. Repository reference files are not installed in Python artifacts.

Repository Files
----------------

The repository keeps these reference assets:

- ``data/hg38.chrom.sizes``
- ``data/hg38_annotated_collapsed_TSS_PAS_1kb.bed``
- ``data/covered_targets_Twist_Methylome_hg38_annotated_collapsed.bed``
- ``data/manifest.json`` and companion ``data/*.npy`` / ``data/*.npz`` arrays
  for model-level power analysis.

The model-power arrays remain in GitHub so the Streamlit deployment can load
``data/manifest.json`` and its companion files from the checkout. They are
excluded from Python wheels and source distributions because of their size.
The 24 tracked array files occupy about 520 MB (496 MiB).

Model-Power Array Layout
------------------------

Keep the manifest and its companion files together in this layout:

.. code-block:: text

   data/
   |-- manifest.json
   |-- cpg_index.npz
   |-- cpg_mean.float32.npy
   `-- std_by_depth/
       |-- depth_5_mean.float32.npy
       |-- depth_5_CI.float32.npy
       `-- ...

The current manifest contract is ``format_version: 1`` and contains 3,771,981
CpG rows. It records the supported depth labels, statistics, dtype, row count,
and relative filenames. Do not move individual arrays or combine files from
different reference releases. Source pickles used to generate the arrays are
not required by the loader or Streamlit app.

Installed API callers must pass ``reference_dir`` to
``analysis.model_power.load_default_model_power_reference`` or set
``CFTK_MODEL_POWER_DATA``. Resolution order is:

1. explicit ``reference_dir``;
2. ``CFTK_MODEL_POWER_DATA``;
3. repository-local ``data/`` when running from a checkout.

The loader validates the manifest and all files required for the requested
depths and SD statistics before reading arrays.

Two Power-Analysis Inputs
-------------------------

The power workflows use different reference inputs and answer different
questions:

``cftk power``
   Uses the user-supplied pickled ``reference_data.cpg_std`` table. It performs
   CpG-level analytical power calculations over sample-size, methylation-effect,
   and depth scenarios.

Model-development API and Streamlit calculator
   Use ``data/manifest.json`` and its NumPy arrays. They simulate the complete
   cross-validated biomarker-discovery pipeline and matched null calibration.

The ``cpg_std`` configuration path is not a substitute for the calculator
array directory, and the calculator arrays are not consumed by ``cftk power``.

The example config still uses placeholder paths. Replace the placeholders with
paths that are valid in your environment.

User-Supplied Files
-------------------

Most full workflows require additional references:

- ``genome_fa``: hg38 FASTA used by alignment, methylation extraction, and some
  QC steps.
- ``genome_2bit``: 2bit genome used by some fragmentomics workflows.
- ``ctcf_bed``: CTCF regions for cleavage.
- ``blacklist`` and ``gap``: excluded regions for DELFI-style features.
- ``bins``: genomic bins for DELFI-style features.
- ``cpg_std``: pickled CpG standard-deviation table used only by the legacy
  ``cftk power`` workflow.

Coordinate Safety
-----------------

Do not mix genome builds or chromosome naming conventions. Keep FASTA,
chromosome sizes, BED files, bins, blacklist, and gap files on the same build
and naming scheme.
