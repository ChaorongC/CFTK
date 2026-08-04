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

**Status:** TODO, starts after managed reference acquisition is stable.

Why this is pending: A diagnostic command needs a stable tool, profile, and
acquisition contract. Implementing it earlier would duplicate evolving
initialization logic and provide incomplete assurances.

Implementation:

- Check the Python package and supported Python version.
- Check required external executables and capture their versions, including
  Trim Galore, bwa-meth, Sambamba, samtools, Picard, MethylDackel, bedtools,
  MultiQC, and workflow-specific optional tools.
- Validate the selected profile manifest, hashes, assay/genome match, target
  BED coordinates, chromosome sizes, bwa-meth indexes, `.fai`, and `.dict`.
- Validate the sample sheet, input readability, output-directory writability,
  available disk space, and config/lock consistency.
- Provide concise human output plus machine-readable JSON and meaningful exit
  codes for schedulers and CI.

Completion criteria:

- `cftk doctor` returns zero only when the selected workflow can start with all
  required inputs, tools, and references.
- Every failed check identifies the exact component and an actionable remedy.
- Tests cover missing, incompatible, corrupt, and permission-denied states.

## 3. Beginner `cftk run`

**Status:** TODO, requires defined failure and resume semantics.

Why this is pending: The current `run-all` command can continue after some
failures and advanced workflows require different references. A beginner
command must not imply reliable completion when only part of a workflow ran.

Implementation:

- Define the beginner default as core processing plus QC; keep differential,
  DMR, fragmentomics, MESA, and power workflows explicit.
- Validate configuration, inputs, tools, and required profile components before
  launching work.
- Define fail-fast behavior, checkpoint validity, resume behavior, dry-run
  output, cancellation, and final success criteria.
- Record resolved paths, config/lock hashes, command lines, external-tool
  versions, and output checksums in run provenance.
- Return nonzero when any required stage fails or its expected artifacts are
  missing. Do not implement this as a silent alias for `run-all`.

Completion criteria:

- A complete run is distinguishable from partial, failed, skipped, and resumed
  runs through both exit status and machine-readable provenance.
- Dry-run and small integration fixtures verify the exact stage graph and
  command construction.
- Interrupted runs resume only from scientifically valid checkpoints.

## 4. Real-Data End-To-End Validation

**Status:** TODO, requires external tools and an approved validation dataset.

Why this is pending: Unit tests and mocked commands validate orchestration but
cannot establish biological equivalence or real-tool compatibility.

Implementation:

- Select an approved small Twist Human Methylome paired-FASTQ dataset and pin
  its inputs, reference profile, expected tool versions, and checksums.
- Execute reference preparation, trimming, bwa-meth alignment, Sambamba
  duplicate marking, Picard metrics, MethylDackel extraction, CpG merge, and QC.
- Compare commands and key outputs with the Twist technical note, including
  read groups, mapping filters, target metrics, `--mergeContext`, minimum depth,
  and OT/OB handling.
- Quantify Sambamba-versus-Picard duplicate-marking agreement rather than
  assuming equivalent output.
- Record runtime, memory, artifact checksums, expected ranges, and any accepted
  version-specific differences.

Completion criteria:

- Every external command exits successfully and produces its expected,
  nonempty artifacts.
- Target coverage, alignment, duplicate, methylation, and QC outputs pass
  predefined scientific and structural checks.
- The complete validation recipe and provenance are reproducible in a clean
  environment before a production release is tagged.
