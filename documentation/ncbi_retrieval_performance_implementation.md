# NCBI Retrieval Performance and Correctness: Implementation Report

Companion to `ncbi_retrieval_performance_plan.md`. That document is the plan;
this one records what was built, what was measured, and where the result
deviates from the projection.

Status: phases P0 through P5 implemented. Two NCBI requests were issued, both
for the P0 invariance test the plan authorizes. No full re-fetch was performed
and the frozen snapshot in `data/01-raw/` is untouched.

Test suite: `python -m unittest discover -s tests -q` — 135 tests, all passing.

---

## 1. Measured results

All figures below are measurements taken against the real artifacts, not
estimates. Where the plan carried an estimate, the measurement is compared to it.

### 1.1 P0 `rettype=gp` invariance test (2 NCBI requests)

Executed by `scripts/ncbi_rettype_invariance_check.py` over the first 100 UIDs of
the frozen snapshot (`1000250755..1011828871`).

| Fetch | Parameters | SHA-256 | Bytes |
| --- | --- | --- | --- |
| A | `retmode=xml`, no `rettype` | `537d492c101433aebc2c7bbad64103cfb35af1ac389da55b52c261e4060de3c3` | 919,109 |
| B | `retmode=xml`, `rettype=gp` | `537d492c101433aebc2c7bbad64103cfb35af1ac389da55b52c261e4060de3c3` | 919,109 |

**Primary criterion `sha256(A) == sha256(B)`: PASS.** `rettype=gp` is a no-op
today, so pinning the documented value costs nothing and removes the dependency
on undocumented default resolution. Acceptance criterion 7 is met.

**Secondary signal.** Both payloads differ from the frozen
`batches[0].xml_payload_sha256` (`c50dc7d3...`). Since A and B agree with each
other, this is record revision upstream since April 2026, not a format change.
Recorded, not treated as a failure — and it is direct evidence for open
decision 3: the NCBI records behind this snapshot have moved.

**Unplanned confirmation.** The first fetch hit two consecutive HTTP 502
responses, retried, and recovered. The failure controls were exercised against
live NCBI rather than only against mocks.

### 1.2 P1 consolidation is byte-identical at full scale

The consolidated snapshot hash is part of the artifact contract, so the
streaming writer must reproduce ElementTree's serialization exactly. Feeding the
real 363 MB artifact's 41,345 records back through the streaming writer:

| Metric | Value |
| --- | --- |
| Records | 41,345 |
| Bytes written | 363,034,184 (identical to source) |
| Streamed SHA-256 | `ada0563932d20f54f68c614df3c40d4fa3520038ad438a41cc471a1eacb970bf` |
| Frozen `xml_file_sha256` | `ada0563932d20f54f68c614df3c40d4fa3520038ad438a41cc471a1eacb970bf` |
| **Identical** | **Yes** |
| Peak RSS above baseline | 3.5 MB |
| Elapsed | 28.5 s |

The unit suite additionally pins byte-equality against the original in-memory
algorithm across entity escaping, non-ASCII text, self-closing elements, root
attributes, and record tail whitespace.

### 1.3 P1 memory: measured, not estimated

The plan estimated the consolidation peak at 2.5-3.5 GB. The dominant term in
that estimate is one full ElementTree of the consolidated document, which is
what the previous validation built — twice per creation and twice per reload.

| Path | Peak RSS above baseline | Elapsed |
| --- | --- | --- |
| Previous whole-document validation, **one** pass | **1.79 GB** | 10.75 s |
| New streaming validation, full reuse path including the manifest hash check | **9.4 MB** | 11.35 s |

A ~190x reduction in peak memory at equal wall time for a single pass, and the
new path performs one pass where the old performed two. The plan's 2.5-3.5 GB
creation estimate is confirmed as well founded: creation additionally held the
raw batch bytes (363 MB) and the `tostring` output (363 MB) alongside the tree.

### 1.4 Observed EFetch latency

From the invariance test, `batch_size=100`, GBSeq XML:

| Request | Round trip | Bytes |
| --- | --- | --- |
| A (after two 502 retries) | 6.52 s | 919,109 |
| B | 4.59 s | 919,109 |

Two samples are not a latency distribution, but they are the first direct
measurement of `t_gen` for this workload and they are consistent with the ~30
minute baseline for 414 sequential batches. They also support the plan's model
conclusion in section 2: at roughly 46 ms of server-side rendering per record,
raising `batch_size` removes per-request overhead but not rendering time, so P3
alone should land near the projected 9%, and concurrency is the only first-order
lever on the remote side. Phase P0 telemetry now records the full distribution
on any real run.

---

## 2. What was built, by phase

### P0. Instrumentation and safety prerequisites

New module `src/pago_pipeline/ncbi_telemetry.py`.

- `NCBIRetrievalTelemetry` collects, per stage (`uid_initial_search`,
  `uid_pages`, `xml_batches`): request count, retry count, reused-batch count,
  response bytes, cumulative regulatory sleep, stage wall time, and latency
  summaries (median, p95, max) split into **connect / time-to-first-byte /
  response-read**. Failure counters are categorized as `rate_limit`,
  `http_429`, `http_5xx`, `timeout`, `deadline`, `response_validation`, `other`.
- Local phases are timed separately: `consolidated_xml_write`,
  `consolidated_xml_validation`, `latest_directory_copy`.
- Process metrics come from `psutil` when available: observed maximum RSS and
  disk read/write deltas. Telemetry never raises; instrumentation must not be
  able to fail the run it observes.
- The summary is written into the XML manifest under `retrieval_telemetry` and
  returned in the resolve payload. The manifest copy necessarily predates the
  `latest/` copy phase, because the manifest is written first; the payload copy
  covers the whole run.

**Connect versus TTFB.** `urllib` exposes no connection-establishment hook, so
`connect_seconds` is populated only by the keep-alive transport (section P5).
TTFB and read time are always split, which is what distinguishes "NCBI is
rendering" from "the link is slow".

**F3 — rate-limit responses are transient.** `_payload_indicates_ncbi_rate_limit`
matches the documented `{"error":"API rate limit exceeded","count":"11"}` body
both as raw text and as parsed JSON, and raises `NCBIRateLimitError`, a subclass
of `NCBITransientFetchError`. The retry backoff for that error is floored at
`rate_limit_backoff_seconds` (default 5.0), because backing off by the ordinary
transient delay would re-enter the same violation. This was the hard gate on P5.

**F10 — bounded batch context.** `_build_xml_batch_context` now emits
`protein_uids=<first>..<last> (n=N)`. The full UID list, request URL, and
response headers move behind `verbose_batch_logging`.

**Live progress and projected finish.** Per-request latency alone leaves the
reader to do arithmetic before deciding whether to let a long run finish, and
under concurrency the per-request lines from several workers interleave and
arrive out of order, which buries the answer entirely. Those lines moved behind
`verbose_batch_logging`, and one self-rewriting line now carries the run:

```
[██████░░░░░░░░░░░░░░░░░░]  76/612  12.4% | #78 23.6s 0.97 MB | 73.6 MB at 296 KB/s | 4m15s elapsed | ~29m57s left | ends 19:09:10
```

Both retrieval stages share one reporter, so a reader does not have to learn
two formats. They differ only in the unit they advance by: the XML stage counts
batches, the UID stage counts identifiers, because that is the quantity the
reader is actually waiting on rather than an artifact of how it was paged.

```
[███████████████░░░░░░░░░] 30,000/46,400  64.7% | page 3 0.3s 0.08 MB | 0.2 MB at 146 KB/s | 1.6s elapsed | ~0.9s left | ends 19:23:09
```

Design notes worth keeping:

- The projection uses completion throughput rather than per-request latency,
  because the two diverge as soon as requests overlap.
- Batches reused from the workspace count as progress but are excluded from the
  rate. They cost no request, so including them would project a finish time no
  real request could meet.
- Rate and projection are withheld until at least one second has elapsed.
  Dividing by a near-zero interval produced confidently wrong numbers such as
  `26,806,283 KB/s`.
- The line is rewritten in place only when the stream can support it, detected
  as a TTY or a live `ipykernel`. A redirected stream gets one line per batch
  instead of accumulating carriage returns.
- Bar glyphs fall back from `█░` to `#-` when the output encoding cannot
  represent them.
- Retry and failure messages are never suppressed: they clear the progress
  line, print on their own row, and the bar is redrawn beneath them.
- Retry count appears on the bar as soon as it is non-zero.

The stage ends with total elapsed time and aggregate transfer rate.

### P1. Streaming persistence and validation (0 NCBI requests)

New module `src/pago_pipeline/ncbi_xml_stream.py`.

- **Spill on arrival.** Each validated batch payload is written atomically to
  the workspace as soon as it is validated; only its SHA-256, counts, and file
  path stay in memory. `NCBIProteinXmlBatchFetchResult` carries either
  `xml_payload_bytes` or `xml_payload_file_path`;
  `read_ncbi_xml_batch_payload_bytes` reads whichever is present.
- **One validation pass per batch.** `validate_xml_batch_payload` checks
  well-formedness, record tag, record count, and the **exact protein UID
  multiset** in a single pass. This is stronger than the previous count-only
  check and is what makes an omitted or substituted record diagnosable at the
  batch that produced it. The failing UIDs are named in the error.
- **Streaming consolidation.** `write_consolidated_xml_document` walks each
  batch with `iterparse` and re-serializes each record straight into the
  destination file. Records are emitted one event behind their end tag, because
  ElementTree assigns an element's tail text only after the following event, and
  that tail is part of the byte-exact output. The root open tag is derived by
  serializing an empty element and converting its self-closing suffix, which
  keeps attribute escaping inside ElementTree rather than reimplementing it.
- **Hash during write.** The SHA-256 is accumulated over the bytes as they are
  written, removing the separate hashing pass over the finished file.
- **One validation pass over the published file.**
  `validate_consolidated_xml_file` checks root tag, record tags, record count,
  and the exact UID multiset in a single streaming pass, replacing the previous
  separate structure pass, UID pass, and hash pass.
- **No reload after publication.** `save_ncbi_protein_xml_snapshot` returns a
  `SavedXmlSnapshot` carrying the manifest, paths, validated record count, and
  hash. `resolve_ncbi_protein_xml_snapshot` builds its return payload from that
  instead of re-reading and re-validating what it just wrote.
- **Notebook 01 stops rehashing.** The consolidated hash is read from the
  manifest, with `VERIFY_XML_FILE_SHA256_BY_REHASHING` available for an explicit
  independent check. Cell 10 states which mode produced the value rather than
  printing a tautology.

**Pass count, stated honestly.** Passes over the *consolidated file* drop from
six reads plus a copy to one read plus a copy. Passes over *batch* data remain
two — the fetch-time validation and the `iterparse` consolidation that P1.3
prescribes — so the total is three full parses, not the one-to-two the summary
table projected. Dropping to two would require abandoning either per-batch
validation, which F6 depends on, or `iterparse` consolidation, which P1.3
requires. The memory result is unaffected.

**The 363 MB copy into `latest/` is retained.** It is listed in F5 as an
observation but is not among P1's numbered items, and removing it would mean
either abandoning the `latest/` artifact contract or hard-linking, which would
let a write to `latest/` mutate the immutable snapshot.

### P2. UID path via History (F1, F11)

`fetch_ncbi_protein_uid_snapshot` now executes the query **once** with
`usehistory=y` and reads every page back from that frozen History set with
EFetch `rettype=uilist`, `retmode=text`, `retstart`, `retmax`.

- Requests for the reference dataset: **43 → 6** (-86.05%), regulatory sleep
  5.46 s → 0.65 s.
- The primary justification is snapshot isolation, not time. Re-running the
  query per page is not atomic: an index change between pages can duplicate or
  silently drop UIDs, and `deduplicate_uids=True` masks duplicates while
  omissions are undetectable. It also stops leaving 42 orphaned History entries
  per run.
- Ordering: History order differs from ESearch relevance order, but the snapshot
  normalizes with `sort_uids=True`, so `protein_uids_sha256` is unchanged for an
  unchanged result set.
- **F11**: the UID path now shares the XML path's failure controls. Both stages
  run through one `_execute_ncbi_request_with_controls` loop, so the two cannot
  drift apart: network timeout, absolute request deadline, bounded retries with
  exponential backoff, transient/permanent classification, rate-limit handling,
  and telemetry are identical. Page payloads are validated for numeric content
  and for not exceeding the requested page size.
- `page_size` defaults to 10,000 and is rejected above it, the documented EFetch
  paging ceiling.

### P3. Larger XML batches

`batch_size` default 100 → **500**, sequential. Requests **414 → 83** (-79.95%),
regulatory sleep 53.82 s → 10.79 s. `batch_size` above 10,000 is rejected.
POST is already selected at ≥200 UIDs, so no transport change was needed. The
EPost variant remains deferred, as the plan directs.

### P4. Selective resume

New module `src/pago_pipeline/ncbi_batch_workspace.py`.

- `batch_plan.json` records the plan identity (`protein_uids_sha256`, UID count,
  `batch_size`, database, `retmode`, `rettype`), the request policy for audit,
  and every planned batch range with its UID hash.
- A batch is reused only when its recorded identity, its stored payload hash,
  and a **full structural and exact-UID revalidation of the payload** all agree
  with the current plan. Anything else is deleted and re-fetched.
- An incompatible plan resets the workspace before any request is issued.
- The workspace lives at `<xml_snapshot_root>/.batch_workspace`, is git-ignored,
  and is purged after a successful publication. A failed run keeps it, which is
  the entire point: the rerun re-fetches only what is missing.
- When no workspace directory is given, an ephemeral spill directory is used and
  removed if the fetch raises.

This does not reduce ideal-run time. It reduces expected time given the
five-retry and circuit-breaker regime, and it is what makes `batch_size=500`
safe to adopt: a failure near batch 80 of 83 no longer discards 80 batches.

### P5. Bounded concurrency (default off)

- `max_concurrent_requests` defaults to **1**. Section 6.4 of the plan requires
  the batch-size change to be measured sequentially before concurrency is
  enabled, so concurrency is opt-in rather than a silent default.
- `NCBIRequestStartRateLimiter` reserves slots keyed on request **start**
  instants under one shared lock, replacing per-worker `time.sleep`. Defaults:
  8 starts/s with an API key, 2.5 without.
- `NCBIAdaptiveConcurrencyGovernor` bounds in-flight requests and, on a
  rate-limit, 429, 5xx, timeout, or deadline signal, halves the request rate and
  drops concurrency one step. Recovery takes one step at a time and only after a
  run of clean responses.
- **Ordered assembly.** Results are written into a slot addressed by
  `batch_index`, never appended. Completion order under concurrency is not plan
  order, and `xml_file_sha256` is defined by plan order alone. A test drives
  batches to complete out of order and asserts both the ordering and that a
  concurrent run publishes the same hash as a sequential one.
- **ContextVar propagation.** `copy_context()` is taken per submission inside
  the `_configured_ncbi_entrez_urlopen` block. A single `Context` cannot be
  entered concurrently, so one copy per task is required, not one shared copy. A
  test asserts worker threads observe the configured SSL context; without this,
  a custom CA would be silently dropped on TLS-intercepting networks.
- A permanent failure in any batch aborts the run and propagates the original
  exception; queued batches short-circuit rather than issuing further requests.

**F7 — connection reuse, opt-in, no new dependency.** `reuse_http_connection`
routes requests through `_NCBIPersistentHttpsTransport`, a per-thread keep-alive
`http.client` connection carrying the same `ssl.SSLContext` as the urllib path.
This resolves open decision 2 without adding `httpx` or `requests`: `requests`
maps `verify=` onto a CA *file* and would not preserve `ssl_ca_directory`
semantics, whereas `http.client` takes the existing context object unchanged. A
stale keep-alive connection reconnects once; an abandoned response body discards
the connection; the transport also supplies the `connect_seconds` measurement
the urllib path cannot.

---

## 2b. Two defects found while running the real workflow

Neither was in the plan. Both were found by running the pipeline on real data.

### Publishing `latest/` could destroy it

`_replace_latest_directory` copied the new artifacts to a temporary directory,
then deleted `latest/`, then renamed the temporary directory into place. On
Windows a directory containing a file another process holds open cannot be
removed, and `shutil.rmtree` deletes what it can before failing. Two real runs
therefore left `latest/` with the payload and protein UID list still present and
`manifest.json` gone, which no reader can validate.

The holder was identified: a PyCharm process indexing the 363 MB
`protein_records.xml`. Opening the file for reading or writing succeeded, but
renaming it failed with `WinError 32` and renaming the directory with
`WinError 5`, which is the signature of a handle opened without
`FILE_SHARE_DELETE`. This is the same failure class as commit `1d71766`, which
fixed it for memory-mapped `.npy` artifacts elsewhere in the pipeline.

The publication order is now: stage the new artifacts, move the existing
`latest/` aside, put the new one in place, then delete what was moved aside. A
failure to move the old directory aside now aborts **before** anything is
destroyed, and a failure to install the new one restores the old. Leftovers from
an interrupted publication are swept by the next one. The immutable snapshot is
never involved, so a publication failure has always been recoverable — the
change is that `latest/` now survives it.

Operationally: exclude `data/` from IDE indexing, or the abort will simply
repeat.

### A truncated response was classified as a permanent failure

A real run at `batch_size=500` aborted on its first batch with
`IncompleteRead(4112725 bytes read)`: NCBI announced a body length, sent 4.1 MB,
and the connection closed early. Retrying is the only correct response, because
the request was accepted and the server had begun answering it.

`_is_transient_ncbi_fetch_error` missed it. `http.client.IncompleteRead` derives
from `HTTPException`, not from `OSError` or `URLError`, and its message contains
none of the wording the text heuristics matched. Two neighbouring cases had the
same hole: `RemoteDisconnected` is a `ConnectionResetError` whose message reads
"Remote end closed connection without response", and `ssl.SSLEOFError` was
claimed by `_is_ncbi_tls_configuration_error` because it subclasses `SSLError`,
so a mid-stream TLS drop was reported as a certificate-authority
misconfiguration.

All three are now transient: `http.client.HTTPException`, `ConnectionError`, and
the two transport-level TLS errors. `SSLCertVerificationError` still classifies
as a configuration problem. Truncation is counted under its own telemetry kind,
`truncated_response`, and counts as an overload signal for the concurrency
governor, so a link that keeps dropping large bodies reduces pressure instead of
retrying into the same wall.

This was latent before P3 and became reachable because of it: a 4.8 MB body
spends roughly nine times as long exposed to a mid-stream drop as the 0.9 MB
body a `batch_size=100` request produced.

---

## 3. Contract and compatibility

- `snapshot_format_version` is bumped to **1.1** for both artifact types.
  `SUPPORTED_XML_SNAPSHOT_FORMAT_VERSIONS = {"1.0", "1.1"}`, so the existing
  frozen snapshot still loads; an unknown version is now rejected explicitly
  rather than being silently misread. This resolves open decision 1.
- Manifest additions are additive: `rettype`, `reused_batch_count`,
  `fetched_batch_count`, `request_policy`, `retrieval_telemetry`, and
  `reused_from_workspace` per batch. UID manifests gain
  `uid_retrieval_strategy`, `esearch_request_count`, `efetch_request_count`,
  `fetch_timeout_seconds`, `request_deadline_seconds`, `retrieval_telemetry`.
- Artifact set, atomic publication, immutable snapshot directories, upstream
  provenance fields, and all existing timeout, deadline, retry,
  partial-response, and circuit-breaker controls are unchanged.
- `NCBIProteinXmlBatchFetchResult` field order changed so the payload can be
  optional. All construction in the repository is keyword-based.

---

## 4. Acceptance criteria

| # | Criterion | Status |
| --- | --- | --- |
| 1 | UIDs from the same History set equal the frozen set | Structurally enforced; needs a real UID run to confirm end to end |
| 2 | Requested UIDs equal UIDs in the returned XML, exactly | Done — enforced per batch and again on the published file |
| 3 | Manifests, hashes, provenance, immutability, atomic publication preserved | Done — verified byte-for-byte against the real 363 MB artifact |
| 4 | Existing failure controls preserved | Done — shared control loop, existing tests unchanged in intent |
| 5 | Zero 429 responses in the full trial | Pending the full trial |
| 6 | Reduction in time, peak memory, and full passes | Memory and passes measured; wall clock pending the full trial |
| 7 | `sha256(A) == sha256(B)` in the rettype invariance test | **Done — PASS** |
| 8 | Explicit prior-manifest compatibility or documented version bump | Done — bumped to 1.1, 1.0 still accepted |
| 9 | `python -m unittest discover -s tests -q` passes | Done — 135 tests |

Criteria 1, 5, and 6's wall-clock component require the trials in section 5.
They are not claimed as met.

---

## 5. Running the trials

Per the plan's measurement protocol: hold one variable per trial, use the
existing immutable UID snapshot for all XML trials, and run the full
confirmation inside NCBI's recommended window for large jobs — weekends, or
21:00-05:00 US Eastern on weekdays.

```bash
python scripts/ncbi_rettype_invariance_check.py
```

For the XML trials, edit CELL 3 of `notebooks/01_ncbi_protein_xml_snapshot.ipynb`
and set `XML_SNAPSHOT_MODE = SnapshotMode.create_new`:

1. **P3 sequential.** `XML_BATCH_SIZE = 500`,
   `XML_MAX_CONCURRENT_REQUESTS = 1`, `XML_REUSE_HTTP_CONNECTION = False`.
   Record `retrieval_telemetry` from CELL 6.
2. **F7.** Same, with `XML_REUSE_HTTP_CONNECTION = True`. The delta is the
   connection-setup cost, and `connect_seconds` becomes populated.
3. **P5.** `XML_MAX_CONCURRENT_REQUESTS = 2`, then 4 only if the previous trial
   shows idle latency and `failure_counts` with zero `http_429` and zero
   `rate_limit`.

Start with a ~5,000 UID trial before the full confirmation run. The workspace at
`<xml_snapshot_root>/.batch_workspace` makes an interrupted trial cheap to
resume.

---

## 6. Open decisions after implementation

| # | Decision | Resolution |
| --- | --- | --- |
| 1 | Manifest compatibility | Resolved: `snapshot_format_version` 1.1, 1.0 still accepted, unknown versions rejected |
| 2 | New HTTP dependency | Resolved without one: stdlib `http.client` keep-alive transport, opt-in, identical CA semantics |
| 3 | Scientific refresh | Still open, and now better informed: the invariance test shows the underlying records have been revised since April 2026 |
| 4 | `GBSeq_references` | Unchanged: 20.6% of payload, no consumer, no API parameter to suppress it |
