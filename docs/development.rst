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
- internal validation artifacts for the default Sambamba duplicate-marking
  implementation. These comparisons are not beginner workflow gates.

The pinned managed default profile, local reference profiles, schema-v2 project
initialization, doctor diagnostics, and manifest-backed beginner workflow are
available now. Legacy configuration and alternative-tool compatibility remain
available only through the expert commands.

Internal Validation
-------------------

The historical cohort audit and Sambamba-versus-Picard comparison are retained
for maintainers and advanced users. They are technical diagnostics, not a
required beginner step and not a claim that duplicate-marking tools are
universally interchangeable.

Every terminal ``cftk run`` attempt now creates that evidence bundle under its
attempt directory. The source-checkout helper
``scripts/validation/summarize_workflow_run.py`` remains as a compatibility
wrapper for historical ``run.json`` files and release records. Keep evidence
tables outside Git; the beginner stage contract in
:doc:`user_guide/beginner_run` remains the authoritative list of required
outputs, while the reporter checks and visualizes that contract.

.. toctree::
   :maxdepth: 1

   user_guide/validation_acceptance
