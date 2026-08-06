# CFTK Roadmap

The local schema-v2 initializer, Twist-aligned processing foundation, and
managed reference acquisition are implemented. Remaining sections must not be
advertised as supported until their completion criteria are met.

## 1. Managed Reference Acquisition

**Status:** COMPLETE. The pinned production profile and verified acquisition
path are implemented.

The production profile uses accession-versioned NCBI/UCSC genome resources and
the maintainer-authorized CFTK covered-target BED pinned to a Git commit.

Implementation:

- [x] Define and version a registry schema containing profile ID, profile version,
  assay, genome build, component URLs, byte sizes, SHA-256 checksums, licenses,
  and source attribution.
- [x] Decide which CFTK-owned assets can be hosted in GitHub and which
  third-party genome assets must remain at authoritative sources.
- [x] Confirm target-BED distribution authorization and reference identity, then
  publish a checksummed Twist-compatible Human Methylome hg38 profile release.
- [x] Download into a staging directory, verify every component, validate BED and
  chromosome compatibility, then atomically publish
  `<reference_root>/<profile>/<version>`.
- [x] Refuse unapproved URL declarations, checksum mismatches, incomplete
  profiles, and accidental replacement of an installed immutable version.
- [x] Add retry, interrupted-download, offline, concurrent-install, and corruption
  tests without weakening the existing local-profile mode.

Completion criteria:

- A fresh environment can run `cftk init` with managed mode and obtain the
  pinned default profile without manually supplying component paths.
- Repeated initialization is idempotent and produces the same profile and lock
  hashes.
- Corrupt, incomplete, incompatible, or unlicensed registry entries fail before
  project configuration is accepted.

## 2. `cftk doctor`

**Status:** IMPLEMENTED; real-environment validation remains in section 4.

Implementation followed the managed-reference work so diagnostics could reuse
the stable tool, profile, lock, and acquisition contracts without downloading
or repairing resources.

Implementation:

- [x] Check the Python package and supported Python version.
- [x] Check required external executables and capture their versions, including
  Trim Galore, bwa-meth, Sambamba, samtools, Picard, MethylDackel, bedtools,
  MultiQC, and workflow-specific optional tools.
- [x] Validate the selected profile manifest, hashes, assay/genome match, target
  BED coordinates, chromosome sizes, bwa-meth indexes, `.fai`, and `.dict`.
- [x] Validate the sample sheet, input readability, output-directory writability,
  available disk space, and config/lock consistency.
- [x] Provide concise human output plus machine-readable JSON and meaningful exit
  codes for schedulers and CI.

Completion criteria:

- `cftk doctor` returns zero only when the selected workflow can start with all
  required inputs, tools, and references.
- Every failed check identifies the exact component and an actionable remedy.
- Tests cover missing, incompatible, corrupt, and permission-denied states.

## 3. Beginner `cftk run`

**Status:** IMPLEMENTED AND STRUCTURALLY VALIDATED; biological acceptance
remains in section 4.

The beginner command is intentionally narrower than `run-all`: schema-v2 core
processing plus QC using the validated default toolchain.

Implementation:

- [x] Define the beginner default as core processing plus QC; keep differential,
  DMR, fragmentomics, MESA, and power workflows explicit.
- [x] Validate configuration, inputs, tools, and required profile components before
  launching work.
- [x] Define fail-fast behavior, checkpoint validity, resume behavior, dry-run
  output, cancellation, and final success criteria.
- [x] Extend the implemented exact-command ledger with resolved config/lock
  hashes, external-tool versions, expected-artifact validation, per-attempt
  mirrors, events, and human/machine-readable summaries. Large output checksums
  remain deferred because hashing BAMs adds substantial duplicate I/O; config,
  lock, and run options are hashed now.
- [x] Return nonzero when any required stage fails or its expected artifacts are
  missing. Do not implement this as a silent alias for `run-all`.

Completion criteria:

- A complete run is distinguishable from partial, failed, skipped, and resumed
  runs through both exit status and machine-readable provenance.
- Dry-run and small integration fixtures verify the exact stage graph and
  command construction.
- Interrupted runs resume only from scientifically valid checkpoints.

Reporting integration:

- [x] Move the workflow evidence collector into the installed package and
  generate artifact, command, resource, and sanitized QC summaries from
  ``cftk run``.
- [x] Link generated evidence figures from ``run-summary.html`` for dry,
  failed, complete, and resumed attempts without rerunning successful
  bioinformatics stages.
- [x] Define and test ``complete_with_reporting_error`` for analysis-complete
  attempts whose evidence generation fails; a later run rebuilds evidence only.

## 4. Real-Data End-To-End Validation

**Status:** In progress; clean structural end-to-end validation is complete,
but biological equivalence and production acceptance remain open.

Completed validation used one control and one sALS sample, each deterministically
subsampled to five million paired reads. The managed Twist/hg38 workflow ran
through trimming, bwa-meth alignment, Sambamba duplicate marking, Picard target
and alignment metrics, MethylDackel extraction, CpG merge, and QC. Slurm job
`54985221` exited zero, produced a 420,435-CpG by two-sample matrix, and passed
the scripted structural assertions. A separate read-only doctor audit covered
21 controls and 19 sALS BAMs and exposed expected legacy-reference, stale-index,
and incomplete-read-group issues.

This establishes real-tool compatibility for the tested path. It does not
establish biological equivalence or cohort-wide output quality. Sambamba is
the default duplicate-marking implementation; Picard comparison is retained
as an internal advanced diagnostic only.

The internal Phase 16 diagnostic provides a reproducible technical comparison:
one preserved control and one preserved sALS 5-million-pair alignment were
processed by both Sambamba and Picard from the same pre-markdup BAM. Structural
checks passed, duplicate classification agreement was 99.994% for both
samples, Twist target metrics matched at reported precision, and CpG overlap
was 100%. This result is retained for maintainers and does not create a
beginner-facing equivalence gate or change the Sambamba default.

The final beginner-run validation used a clean workflow-profile environment
and the same deterministic two-sample inputs. Slurm job `55026503` exited zero
after recovering from a deliberately preserved failed attempt: 13 partial
trimming artifacts were quarantined, all seven stages completed, all 24
external commands exited zero, 56 output/report rows and 16 figure rows were
validated, and all 95 run-summary links resolved. A second `cftk run` resumed
all seven stages and executed zero external commands. The resulting matrix has
420,435 CpGs and two sample columns.

Implementation:

- Preserve the approved smoke inputs, managed reference profile, tool versions,
  checksums, external validation recipe, and generated reports outside Git.
- Repeat the complete workflow from a clean project so the command ledger also
  includes stages that were reused from checkpoints during the successful retry.
- Compare commands and key outputs with the Twist technical note, including
  read groups, mapping filters, target metrics, `--mergeContext`, minimum depth,
  and OT/OB handling.
- Keep the default workflow on Sambamba; retain Picard as an explicit advanced
  `process.duplicate_marking_tool` option and keep any comparison internal.
- Record runtime, memory, artifact checksums, expected ranges, and any accepted
  version-specific differences.

Completion criteria:

- Every external command exits successfully and produces its expected,
  nonempty artifacts.
- Target coverage, alignment, duplicate, methylation, and QC outputs pass
  predefined scientific and structural checks.
- The complete validation recipe and provenance are reproducible in a clean
  environment before a production release is tagged.

## 5. Known Readiness And Beginner-Workflow Gaps

These limitations are intentionally recorded rather than hidden by the doctor
implementation:

- Historical hg38 BAMs may contain a 455-contig sequence dictionary while the
  managed no-alt profile contains 195 contigs. Shared contig lengths are not
  sufficient: downstream processing must use the exact ordered alignment
  reference, and CFTK does not currently offer the historical profile.
- A BAM path does not declare whether alignment, duplicate marking, and
  extraction prerequisites are complete. Doctor reports missing read groups or
  duplicate-marking provenance, but provenance warnings cannot prove what was
  done outside CFTK.
- BAM indexes can exist but be older than their BAM. Doctor treats this as a
  readiness failure because region queries may be invalid.
- Bismark alternatives remain accepted configuration values but have not been
  validated for output parity with the default bwa-meth/MethylDackel path.
- ``run-all`` continues after some failures and has no scientifically defined
  resume contract. It is an expert compatibility command, not a substitute for
  the implemented beginner ``cftk run`` command.
- Managed reference acquisition and bwa-meth index construction require large
  cache and temporary-storage allocations. Shared HPC users may need to set
  ``CFTK_REFERENCE_ROOT`` to a lab-managed location before initialization.
- Python dependencies are split by the tested beginner core, ``analysis``,
  ``fragmentomics``, and ``web`` command graphs. XGBoost is not required; it is
  an optional, non-default MESA classifier and its Linux wheel pulled a large
  NCCL runtime. Keep each extra validated independently as commands and
  numerical libraries evolve.
- ``process.cores`` is now enforced as one total CPU budget across processing
  and QC. The run manifest records per-stage allocation, doctor checks detected
  scheduler capacity, and nested sample/tool parallelism cannot exceed the
  configured budget. Memory remains site- and tool-dependent and must still be
  requested explicitly from the scheduler.

## 6. Conda Distribution

**Status:** TODO; the current two-step source installation is accepted and
documented, while native Conda/Bioconda distribution remains future work.

- [x] Document the supported source-checkout sequence: create and activate the
  environment, then install CFTK with `python -m pip install .`.

- [ ] Build a versioned Conda recipe from the released CFTK source artifact and
  declare the validated Python and external-tool runtime dependencies in the
  recipe.
- [ ] Publish through an appropriate maintained channel, preferably Bioconda
  for the bioinformatics toolchain, following its recipe, test, license, and
  update-bot requirements.
- [ ] Make a fresh user installation require only a Conda/Mamba command, with
  no user-facing pip step or repository checkout.
- [ ] Add recipe tests for `cftk --help`, `cftk init --help`, `cftk doctor`, and
  a packaged synthetic dry run; retain the separate real-data release gate.
- [ ] Decide how optional analysis, fragmentomics, and web functionality maps
  to Conda packages or documented add-on environments without reintroducing
  XGBoost into the default processing installation.
- [ ] Verify package and dependency licenses and define the tagged-release
  requirements for the Conda source artifact.

Completion criteria:

- A beginner can install the released package and validated default toolchain
  with one Conda or Mamba command and no repository checkout.
- The published Conda artifact reproduces the validated default dependency and
  external-tool graph on a clean supported platform.
- Package metadata, licenses, tests, and provenance are sufficient for a
  reproducible tagged release.
