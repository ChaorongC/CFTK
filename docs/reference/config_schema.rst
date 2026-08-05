Configuration Schema
====================

Schema V2
---------

``schema_version``
   Must be integer ``2``.

``project_name``
   Required project identifier.

``output_dir``
   Optional output root, relative to the config file when not absolute.
   Defaults to ``.``.

``samples``
   Required path to the strict TSV sample sheet, relative to the config file
   when not absolute.

``assay`` and ``genome``
   Default to ``twist_human_methylome`` and ``hg38``. Both must exactly match
   the selected profile manifest.

``reference_mode``
   ``local`` resolves an installed profile. ``managed`` acquires a profile from
   the validated registry and is the default for new projects.

``reference_root``
   Root containing profile directories. ``CFTK_REFERENCE_ROOT`` overrides this
   value at runtime for portability.

``reference_profile``
   A profile ID string or an object containing ``id`` and optional ``version``.
   When version is omitted, exactly one installed version must exist.

``process``
   Optional compact object with positive integer ``cores``,
   ``parallel_samples``, and ``min_depth`` values. Defaults are 20, 1, and 10.
   Duplicate marking uses ``sambamba`` by default. Advanced users may set
   ``duplicate_marking_tool`` to ``picard`` or ``samblaster`` when their
   environment or validation protocol requires it.
   The optional ``picard_java_memory`` value uses JVM size notation such as
   ``8g`` or ``4096m`` and defaults to ``8g``.

Sample Sheet Schema
-------------------

The required columns are ``sample``, ``group``, ``role``, ``input_type``,
``r1``, ``r2``, and ``bam``. The first release supports exactly two groups and
requires one ``control`` role and one ``case`` role. FASTQ samples require R1
and R2; BAM samples require a BAM path. Rows and group insertion order are
preserved.

Resolved Contract
-----------------

CFTK resolves schema v2 internally into the legacy nested keys expected by
existing commands: ``comparison``, grouped ``samples``, ``reference_data``,
``process``, and ``analysis``. The profile ``target_bed`` is used for Picard
coverage only; it does not filter methylation calling or downstream analyses.

Legacy Schema
-------------

Configs without ``schema_version: 2`` retain the established requirements:
``project_name``, ``output_dir``, ``comparison``, ``samples``,
``reference_data``, ``process``, and ``analysis``. Required references are
``genome_fa``, ``genome_2bit``, and ``chrom_sizes``. Legacy behavior and field
names are preserved.
