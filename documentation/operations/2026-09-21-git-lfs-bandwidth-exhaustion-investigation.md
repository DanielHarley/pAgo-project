# Git LFS bandwidth exhaustion investigation

## Date

2026-09-21

## Initial state

GitHub reported that the account had consumed 100% of the included 10 GB Git
LFS bandwidth for the current billing cycle.

The repository contains one current large LFS object used by two checked-in XML
paths:

- `data/01-raw/protein_xml_snapshots/latest/protein_records.xml`
- `data/01-raw/protein_xml_snapshots/snapshots/2026-04-09T00-51-02Z__q_891f443d754c/protein_records.xml`

Both pointers reference the same object:

- SHA-256 OID:
  `ada0563932d20f54f68c614df3c40d4fa3520038ad438a41cc471a1eacb970bf`
- size: 363,034,184 bytes

## What was found

`.github/workflows/ci.yml` used `actions/checkout@v4` with `lfs: true` in
both CI jobs. A sampled successful run for PR #38, workflow run
`35238280184`, showed the checkout step executing
`git lfs fetch origin refs/remotes/pull/38/merge`.

The same run contained both the Linux `environment` job and the Windows
`windows-snapshot-reuse` job. Each job starts from a fresh runner and performs
its own checkout.

The workflow is triggered by both `push` and `pull_request`. The combination
of repeated workflow events and two LFS-enabled jobs therefore created a
repeated download path for the same immutable 363 MB object even when the
changed code and tests did not need the raw XML bytes.

For scale, one hydration is about 0.338 GiB. Thirty fresh hydrations are about
10.14 GiB, enough to exhaust a 10 GiB monthly allowance.

## Remediation

The accompanying change makes routine CI keep LFS pointers unmaterialized and
updates `scripts/verify_raw_data.py` to validate an unmaterialized pointer by
comparing its SHA-256 OID with the manifest's `xml_file_sha256`.

Full byte-level verification remains explicit through `git lfs pull`.

## Scope limits

This record establishes a repository-internal mechanism capable of explaining
the observed bandwidth consumption. The available GitHub connector does not
expose account-level Git LFS billing attribution, so this record does not claim
that every byte in the 10 GB total came from CI. External clones or other LFS
downloads may also have contributed.
