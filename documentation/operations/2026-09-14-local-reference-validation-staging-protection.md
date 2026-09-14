# Local reference and validation staging protection

Date: 2026-09-14

## Purpose

Prevent local, not-yet-reviewed reference-tree-construction and
profile/placement-validation files from being accidentally staged into an
unrelated Git change.

Historical roadmap reference: B5/B6 (metadata only — see
`../SCIENTIFIC_STATE.md`; not the name of this operation).

## Git baseline

`master` = `origin/master` = `7f991f1d78d2fb35f7acaeed8e9b799f06f86058` at the
start of this operation.

## Initial state

Fourteen local files existed inside the repository working tree:

```text
documentation/pago_annotation_ontology.md
notebooks/11_reference_profiles_and_tree.ipynb
notebooks/12_profile_validation.ipynb
src/pago_pipeline/clade_hmm_build.py
src/pago_pipeline/clade_hmm_build_snapshot.py
src/pago_pipeline/midpiwi_reference_tree.py
src/pago_pipeline/midpiwi_reference_tree_snapshot.py
src/pago_pipeline/pago_clade_placement_validation.py
src/pago_pipeline/pago_domain_profile_validation.py
tests/reference_builder_fixtures.py
tests/test_midpiwi_reference_tree.py
tests/test_pago_clade_placement_validation.py
tests/test_pago_domain_profile_validation.py
tests/test_reference_hmm_builders.py
```

Git state for all fourteen at the start of the operation:

- untracked
- not ignored
- absent from `master` and from every branch's recorded history, including
  `(feat)-siepe-ready-project` (`git log --all` returned nothing for each path)

A broad `git add -A` could therefore have staged all fourteen into an unrelated
commit.

## Backup audit

A previously known backup location, a Claude Code session scratchpad directory
outside the repository, was inspected before relying on it.

The directory contained:

- `MANIFEST.tsv` with path, size, and SHA-256 for 13 of the 14 files
- three subdirectories mirroring the source layout (`notebooks/`, `src/`,
  `tests/`), all empty

It did not contain any actual copy of file content. It was therefore classified
as an integrity manifest, not a recoverable backup.

13/14 files were represented in that old manifest.
`documentation/pago_annotation_ontology.md` was not represented in it.

## Integrity evidence

For the 13 files represented in the old manifest, the following three-way
equality was confirmed for every one of them:

```text
old_manifest_sha256 == current_local_sha256 == new_backup_sha256
```

For `documentation/pago_annotation_ontology.md`, which was not covered by the
old manifest:

```text
current_local_sha256 == new_backup_sha256
```

No historical continuity before this operation is claimed for that file — only
its state at the time of this operation is recorded.

## Durable backup

A real, byte-for-byte backup was created outside the repository and outside any
temporary/scratchpad storage, preserving the source-relative directory
structure.

Backup root:

```text
%USERPROFILE%\Documents\pAgo-project-backups\
2026-09-14-b5b6-pre-reintegration\
```

That directory has its own `MANIFEST.tsv` (`relative_path`, `size_bytes`,
`sha256_original`, `sha256_backup`) and `README.md` with the copy timestamp in
UTC, source repository path, `master` baseline, and explicit statements that
the backup is not validation and is not a substitute for reviewed repository
history.

Verification performed:

- 14/14 file sizes match between original and backup copy
- 14/14 SHA-256 hashes match between original and backup copy

The backup is not part of Git and is not referenced from any tracked file's
content.

## Repository protection

Fourteen exact-path `.gitignore` entries were added, one per file, under a
single explanatory comment block.

No glob rule was used. An exact path cannot match any file other than the one it
names, while a broad glob could silently capture unrelated future files.

The rules:

- prevent ordinary broad `git add` operations from staging these paths
- do not delete, move, or modify the files
- do not review or validate the files scientifically
- do not themselves constitute a backup; the external copy above is the backup

Being ignored, being backed up, being reviewed, and being scientifically valid
are different facts. This operation established staging protection and a
separately verified backup only.

## Git provenance

- Branch: `(chore)-protect-local-reference-validation-staging-files`
- Commit: `7c9d2433f724c5c772904c40a611afa20212062f`
- Commit title: `(chore) protect local reference and validation staging files from accidental staging`

## Validation

- `git check-ignore -v` resolves all 14/14 paths to the new `.gitignore` rules
- `git add -A -n` does not select any of the 14 paths
- the pre-existing `data/01-raw/**` default-deny rule from PR #29 is unchanged
- `.gitattributes` is unchanged
- `git diff --check` passes
- original-vs-backup hashes remain identical
- `python -m unittest discover -s tests -q` still collects and attempts the
  same local untracked test files as before this change: 106 tests total with
  the same two pre-existing `ModuleNotFoundError: tests.apaz_hmm_test_support`
  import errors. `.gitignore` is a Git mechanism and does not control Python's
  filesystem-based test discovery.

## Scientific impact

None. This operation changes repository staging safety only. It does not
establish, evaluate, or imply the scientific validity of any reference, model,
tree, profile, classifier, or validation procedure represented by the fourteen
protected files.

## Related records

- `../PROJECT_LEDGER.md` — entry dated 2026-09-14, "Protect local reference and
  validation files from accidental staging"
- PR #31 — merged 2026-09-14, merge commit
  `06c52140dfc0b2ce24240df76f7f9acd7cf2ed8c`
- Historical roadmap reference: B5/B6
