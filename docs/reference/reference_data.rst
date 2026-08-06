Reference Data
==============

CFTK needs both repository-hosted region resources and user-supplied genome
resources. Repository reference files are not installed in Python artifacts.

Versioned Processing Profiles
-----------------------------

Schema-v2 projects expose one ``reference_root`` rather than individual paths.
Each local profile has this layout:

.. code-block:: text

   <reference_root>/
   `-- twist_human_methylome_hg38/
       `-- 1.0.0/
           |-- manifest.json
           |-- genome/
           |   |-- hg38.fa
           |   |-- hg38.2bit
           |   `-- hg38.chrom.sizes
           `-- assay/
               `-- covered_targets.bed

The manifest is illustrated by
``examples/reference-profile-manifest.json``. Required components are
``genome_fa``, ``genome_2bit``, ``chrom_sizes``, and ``target_bed``. Additional
workflow-specific components such as ``tss_pas_bed``, ``ctcf_bed``,
``blacklist``, ``gap``, ``bins``, and ``cpg_std`` may be declared in the same
``components`` object and are exposed to downstream commands when present.

Component paths must remain inside the version directory. Optional manifest
``sha256`` values are checked during initialization; CFTK computes every
component hash for ``cftk.lock.json`` even when the local manifest omits it.
Target BED intervals must be valid 0-based half-open BED coordinates on contigs
present in ``chrom_sizes`` and must not exceed chromosome lengths.

Root resolution is ``CFTK_REFERENCE_ROOT``, the project JSON hint, then
``~/.cache/cftk/references``. The environment override is intended for moving
a project between machines.

Managed Registry And Acquisition
--------------------------------

The packaged managed registry has schema version 1. Each profile is keyed by
profile ID and version and records its assay, genome build, and all required
components. Every component must declare:

- a safe relative installed path, byte size, and SHA-256 checksum;
- one or more HTTPS artifact URLs, artifact size and SHA-256 checksum,
  compression type, and an explicit immutability assertion;
- a license name and HTTPS URL; and
- source attribution with a name and HTTPS URL.

Downloads are written under a temporary directory on the reference-root
filesystem. CFTK verifies the downloaded artifact, decompresses it when needed,
verifies the installed file independently, checks target BED/chromosome
compatibility, and atomically publishes the completed version. A stable install
lock serializes concurrent processes. Repeated installation is idempotent;
CFTK refuses to overwrite an installed version that is corrupt or differs from
the current registry entry.

The packaged registry publishes
``twist_human_methylome_hg38`` version ``1.0.0``. It contains:

- the NCBI ``GCA_000001405.15`` GRCh38 no-alt analysis-set FASTA with UCSC IDs;
- a two-column chromosome-size file deterministically projected from NCBI's
  companion FASTA index;
- UCSC's sequence-matched ``hg38.analysisSet.2bit``; and
- the Twist Human Methylome covered-target BED from CFTK commit ``3cea475``,
  distributed with maintainer authorization.

The NCBI FASTA, NCBI FASTA index, and UCSC 2bit resolve to the same ordered 195
contigs. The BED raw URL is pinned to the Git commit and its bytes match the
repository asset. CFTK records and verifies separate transport and installed
SHA-256 values where gzip decompression or FASTA-index projection changes the
byte stream. The profile is Twist-compatible; CFTK does not claim that its FASTA
filename is a byte-identical copy of every vendor resource bundle.

For controlled development, ``CFTK_REFERENCE_REGISTRY`` may point to an approved
registry JSON file. The same schema, HTTPS, size, checksum, license,
source, staging, and compatibility checks apply. This override is not a way to
bypass data-use requirements.

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

Schema-v2 processing resolves the covered-target BED from its profile for
Picard ``CollectHsMetrics``. Legacy source checkouts retain the bundled
fallback. ``process --target-bed PATH`` is an expert one-run override, and
``--skip-picard-metrics`` explicitly disables these metrics.

For the default Twist targeted assay, the same ``target_bed`` also scopes WPS,
occupancy, and DELFI to panel-overlapping reads and intervals in ``auto`` mode.
Those outputs are panel-restricted rather than genome-wide; use the explicit
``--fragmentomics-scope genome`` override only with validated whole-genome
inputs.

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

Running ``cftk init`` prepares bwa-meth converted indexes, ``genome_fa.fai``,
and a conventional sequence dictionary beside the FASTA (for example,
``hg38.dict`` for ``hg38.fa``). Existing complete files are reused. CFTK
requires the generated or existing ``.fai`` contigs and lengths to exactly
match the profile chromosome sizes.

Coordinate Safety
-----------------

Do not mix genome builds or chromosome naming conventions. Keep FASTA,
chromosome sizes, BED files, bins, blacklist, and gap files on the same build
and naming scheme.
