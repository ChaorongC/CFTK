Development
===========

CFTK is under active development. The current package version is 1.0.0 and the
source repository is https://github.com/ChaorongC/CFTK.

Planned Release Work
--------------------

The repository ``TODO.md`` is the authoritative roadmap. Managed references,
``cftk doctor``, and the fail-fast beginner ``cftk run`` are implemented. The
remaining release gates are:

- end-to-end validation with approved Twist data and the real external
  bioinformatics toolchain before a production release is tagged; and
- predefined biological acceptance criteria, including quantified
  Sambamba-versus-Picard duplicate-marking agreement.

The pinned managed default profile, local reference profiles, schema-v2 project
initialization, doctor diagnostics, and manifest-backed beginner workflow are
available now. Legacy configuration and alternative-tool compatibility remain
available only through the expert commands.
