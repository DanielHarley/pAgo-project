# NCBI Retrieval Performance and Correctness Plan

Scope: notebooks `00_ncbi_protein_ids_query` and `01_ncbi_protein_xml_snapshot`,
backed by `src/pago_pipeline/ncbi_api.py` and `src/pago_pipeline/ncbi_snapshot.py`.

Goal: reduce the wall-clock time of NCBI UID and XML retrieval without weakening
the snapshot, manifest, provenance, and validation contracts.

Status: planning only. No source file has been modified. No NCBI request has been
issued during the preparation of this document. All baseline figures below were
derived from the materialized artifacts already present in `data/01-raw/`.

---

## 1. Measured baseline

These are measurements, not estimates.

### 1.1 Snapshot identity

| Item | Value |
| --- | --- |
| Source UID snapshot | `2026-04-06T20-42-12Z__q_891f443d754c` |
| XML snapshot | `2026-04-09T00-51-02Z__q_891f443d754c` |
| Search query | `PIWI[All Fields] AND Bacteria[Organism]` |
| `protein_uids_sha256` | `6772f5e988d76a3756c95f8781e812bcdb1b0d0e465606a76251286beecc89f5` |
| `xml_file_sha256` | `ada0563932d20f54f68c614df3c40d4fa3520038ad438a41cc471a1eacb970bf` |
| `batches[0].xml_payload_sha256` | `c50dc7d313f9a45befdacb9fa9840ace1579b80164b0c4a163a53c50522b0b3d` |

### 1.2 Volume and request counts

| Metric | Value |
| --- | --- |
| Normalized UIDs | 41,345 |
| Consolidated XML | 363,034,184 bytes |
| ESearch requests (`page_size=1000`) | 43 (1 initial + 42 pages) |
| EFetch requests (`batch_size=100`) | 414 |
| Regulatory sleep, UID stage | 42 x 0.13 s = 5.46 s |
| Regulatory sleep, XML stage | 414 x 0.13 s = 53.82 s |

Sleep uses `request_delay_seconds = 0.10` (API key present, per manifest) plus a
mean jitter of 0.03 s from `random.uniform(0.01, 0.05)`. Sleep excludes network
latency, transfer, parsing, retries, and persistence.

### 1.3 Record size distribution

| Statistic | Bytes per record |
| --- | --- |
| Mean | 8,781 |
| Median | 8,109 |
| p90 | 12,237 |
| p99 | 16,874 |
| Max | 116,124 |

The distribution is tight. There is no heavy tail that would make large batches
unpredictable in size.

### 1.4 Batch payload size as a function of `batch_size`

Computed over the materialized XML in file order.

| `batch_size` | Batches | Mean | p90 | Max |
| --- | --- | --- | --- | --- |
| 100 (current) | 414 | 0.88 MB | 1.05 MB | 1.86 MB |
| 250 | 166 | 2.19 MB | 2.57 MB | 3.23 MB |
| **500** | **83** | **4.37 MB** | **5.16 MB** | **5.46 MB** |
| 1000 | 42 | 8.64 MB | 10.13 MB | 10.52 MB |
| 5000 | 9 | 40.34 MB | 50.44 MB | 50.44 MB |

### 1.5 Payload composition

| Block | Size | Share | Consumed downstream |
| --- | --- | --- | --- |
| `GBSeq_feature-table` | 195.1 MB | 53.7% | Yes: `pago_qc.py` reads `feature__region__qual__{region_name,note,db_xref}`, `feature__site__qual__*`, `feature__protein__qual__*` |
| `GBSeq_references` | 74.8 MB | 20.6% | No consumer identified |
| `GBSeq_sequence` | 22.7 MB | 6.2% | Yes: `metadata_to_fasta.py` |
| Remainder | 70.5 MB | 19.5% | Yes |

**Conclusion: payload reduction is not an available lever.** The dominant block is
required, and EFetch exposes no parameter to suppress `GBSeq_references`.
Switching to `rettype=fasta` or ESummary would break the QC evidence layer.

### 1.6 Operating regime

363 MB retrieved over an approximately 30-minute stage implies an effective
throughput near 200 KB/s, roughly 1-2% of any reasonable link. **The stage is
latency-bound, not bandwidth-bound.** Server-side GenBank XML rendering dominates.

---

## 2. Cost model

```
T_xml(b, P) = ceil(N/b) * [ t_conn + t_gen(b) + t_xfer(b) + t_cpu(b) + t_sleep ] / P
```

with `t_gen(b) ~ alpha * b` (server-side rendering, linear in record count) and
`t_conn` constant per request.

Two consequences that determine priority:

1. **Increasing `b` only reduces the `t_conn + t_sleep` term.** The ceiling on that
   gain is `ceil(N/b) * t_conn`, roughly 166 s within a ~30 min stage, or ~9%.
2. **Increasing `P` divides the entire bracket, including `t_gen`.** This is the
   only first-order lever on the remote side.

`alpha` and `t_conn` are unmeasured. Phase P0 exists to fill them in.

---

## 3. Findings

Severity reflects combined impact on correctness, reproducibility, and time.

| ID | Severity | Area | Finding |
| --- | --- | --- | --- |
| F1 | High | Correctness | UID pagination discards the History it obtained |
| F2 | High | Contract | `retmode=xml` without `rettype` is undocumented for `db=protein` |
| F3 | High | Robustness | Rate-limit response is classified as a permanent failure |
| F4 | High | Memory | Consolidation peak estimated at 2.5-3.5 GB |
| F5 | Medium | Local I/O | 7 full parses/reads plus a 363 MB file copy per creation |
| F6 | Medium | Robustness | One suppressed UID aborts the run; no resume |
| F7 | Medium | Throughput | No HTTP keep-alive; ~166 s of TCP+TLS setup |
| F8 | Medium | Throughput | `batch_size=100` yields 414 requests |
| F9 | Medium | Throughput | Serial at ~0.25 req/s against a 10 req/s allowance |
| F10 | Low | Observability | Per-batch context embeds the full UID list |
| F11 | Low | Consistency | UID path lacks the failure controls the XML path has |
| F12 | - | Correction | `latest/` artifacts are intact; a prior claim to the contrary is withdrawn |

### F1. UID pagination discards the History

`fetch_ncbi_protein_uid_snapshot` captures `WebEnv` and `query_key`, prints them,
stores them in the manifest, and never uses them. Pagination calls:

```python
Entrez.esearch(db=..., term=query, usehistory="y",
               retmax=page_size, retstart=page_start_index)
```

Three consequences:

1. NCBI re-executes the full query on each of the 42 pages.
2. Each call creates a new History entry, leaving 42 orphaned WebEnvs per run.
3. **No snapshot isolation.** `retstart` pagination over a re-executed query is not
   atomic. If the Entrez index changes mid-run, UIDs can be duplicated or silently
   skipped. `deduplicate_uids=True` masks duplicates; omissions are undetectable.

Item 3 contradicts the reproducibility premise stated in the function's own
docstring. This is a correctness defect that happens to also cost time.

### F2. Undocumented EFetch format contract

`_open_ncbi_entrez_efetch_once` sends `retmode="xml"` with no `rettype`. Per
NBK25499 Table 1, for `db = nuccore, protein or popset` the defined combinations are:

| Record type | `rettype` | `retmode` |
| --- | --- | --- |
| text ASN.1 | null | text (default) |
| binary ASN.1 | null | asn.1 |
| Full record in XML | `native` | xml |
| TinySeq XML | `fasta` | xml |
| **GBSeq XML** | **`gp`** | **xml** |
| INSDSeq XML | `gpc` | xml |

`rettype=null` with `retmode=xml` is **not in the table**. NCBI currently resolves it
to GBSeq XML, confirmed by the materialized `<GBSet><GBSeq>` payload, but that is
undocumented behavior.

The entire downstream chain (`ncbi_metadata_csv.py` -> `pago_qc.py` -> FASTA ->
SWeeP) depends on the `GBSeq_*` schema. If NCBI changed the default to `native` or
`fasta`, metadata extraction would break. Setting `rettype=gp` is an active
reproducibility debt, not cosmetic.

### F3. Rate-limit response kills the run

NCBI answers rate-limit violations with `{"error":"API rate limit exceeded","count":"11"}`.
If delivered with HTTP 200, the current path is:

```
_xml_payload_contains_server_error_marker() -> False   (no <ERROR> tag)
ET.fromstring() -> ParseError
_is_transient_ncbi_fetch_error(ParseError) -> False
-> RuntimeError("Permanent NCBI XML batch failure; not retrying")
```

Harmless today at 0.25 req/s. **Phase P5 introduces exactly the condition that
triggers it.** This is a hard prerequisite for concurrency.

### F4. Consolidation memory peak

`_build_consolidated_xml_payload` parses every batch into an ElementTree, reparents
all 41,345 records under one root, then materializes the whole document via
`ET.tostring`. Simultaneously live at that moment:

| Component | Estimate | Why retained |
| --- | --- | --- |
| `fetch_result.xml_batches` raw bytes | 363 MB | `fetch_result` is still needed for the manifest |
| Consolidated ElementTree | 1.5-3 GB | ~41,345 records x ~50 elements = ~2M `Element` objects |
| `consolidated_xml_payload` bytes | 363 MB | output of `tostring` |

Estimated peak: **2.5-3.5 GB**. On a typical Windows workstation this crosses the
paging threshold, at which point consolidation becomes swap-bound and its cost
becomes highly non-linear. This is plausibly the dominant local cost and is not
addressed by any remote optimization.

### F5. Repeated full passes over the XML

Per creation:

| # | Pass | Location |
| --- | --- | --- |
| 1 | `_count_xml_batch_records` per batch | `ncbi_api.py` |
| 2 | `_build_consolidated_xml_payload` re-parses every batch | `ncbi_snapshot.py:91` |
| 3 | `_validate_saved_consolidated_xml_snapshot` on the written file | `ncbi_snapshot.py:613` |
| 4 | `_validate_xml_record_uids` on the written file | `ncbi_snapshot.py:617` |
| 5 | `sha256_of_file` | `ncbi_snapshot.py:619` |
| 6 | `shutil.copy2` to `latest/` (363 MB read + 363 MB write) | `ncbi_snapshot.py:645` |
| 7 | `_validate_saved_consolidated_xml_snapshot` on reload | `ncbi_snapshot.py:772` |
| 8 | `_validate_xml_record_uids` on reload | `ncbi_snapshot.py:776` |
| 9 | SHA-256 recomputed by notebook 01 | `notebooks/01_*.ipynb` |

Seven full parses/reads, one full-file copy, one redundant hash.

### F6. Brittle partial-response policy, no resume

`record_count != len(batch_uids)` is classified permanent and not retried. NCBI
silently omits suppressed or removed UIDs. A single dead UID therefore aborts the
entire run, and there is no checkpoint, so all completed batches are discarded.
At `batch_size=500` the blast radius grows fivefold.

### F7. No connection reuse

`Entrez.urlopen` delegates to `urllib.request.urlopen`, which does not keep
connections alive. Each of the 414 requests pays a TCP handshake plus a TLS
handshake (~2 RTT). At ~120 ms RTT this is ~0.4 s/request, roughly **166 s** of pure
connection setup.

### F12. Correction: `latest/` artifacts are intact

A prior assessment claimed the five `latest/` files had been removed and that a
rerun would therefore enter the live-creation path. That is incorrect.

`git status` reports ` M` (modified), not ` D` (deleted), and `git diff` shows
25 deletions and 25 insertions with identical content on a 25-line manifest, the
signature of a CRLF/LF normalization difference. On-disk sizes match the manifest
exactly: `protein_records.xml` is 363,034,184 bytes and both `protein_uids.txt`
files are 493,124 bytes.

**A rerun under `reuse_latest_or_create` will reuse the frozen snapshot.** Choosing
between scientific refresh and reuse remains a legitimate decision, but it is not
forced by the current tree state, and no full re-fetch is required for this work.

---

## 4. Phased plan

Ordering principle: remove correctness and safety defects first; then extract local
gains that cost zero NCBI requests; then reduce remote request count; then add
concurrency last, because concurrency amplifies every unfixed failure mode.

### P0. Instrumentation and safety prerequisites

No performance change. Establishes the measurements the model needs and removes the
blocker for P5.

1. Instrument per stage: initial search, UID pages, HTTP connect, time to first byte,
   response read, regulatory sleep, retries, XML validation, consolidation, write,
   copy to `latest/`, final load.
2. Record: request count, bytes, median and p95 latency, cumulative sleep, 429 and
   5xx counts, peak RSS, disk read/write volume.
3. Split `round_trip_latency_seconds` into connect / TTFB / read. TTFB versus read
   time distinguishes "NCBI is rendering" from "the network is slow", which is what
   determines whether P5 scales linearly or saturates at P=2.
4. **Fix F3**: detect the JSON rate-limit body explicitly and classify it as transient
   with extended backoff.
5. **Fix F2 with an invariance test** (see 4.1 below).
6. Fix F10: truncate `_build_xml_batch_context` to `first..last (n=N)`; put full
   header and UID dumps behind a verbosity flag.

### 4.1 The `rettype=gp` invariance test

Design matters here. Comparing a fresh payload against the April 2026 hash conflates
two hypotheses: a format change, and ordinary record revision since the snapshot was
taken. The discriminating test is a paired fetch:

1. Take the first 100 UIDs of `protein_uids.txt` (`1000250755` first).
2. Fetch A: `retmode=xml`, no `rettype` (current behavior).
3. Fetch B: `retmode=xml`, `rettype=gp`.
4. **Primary criterion**: `sha256(A) == sha256(B)`. Equality proves `rettype=gp` is a
   no-op today and can be adopted without touching the consolidated hash.
5. **Secondary signal**: compare both against `c50dc7d313f9a45befdacb9fa9840ace1579b80164b0c4a163a53c50522b0b3d`.
   A mismatch here with A == B indicates record revision, not a format change, and
   should be recorded rather than treated as a failure.

Cost: 2 requests.

### P1. Streaming persistence and validation (offline)

This phase issues **zero NCBI requests**. It can be developed and tested entirely
against the existing immutable snapshot, which makes it the cheapest phase to
validate and the safest to do early.

1. Persist each batch to an atomic temporary file immediately after transfer, keeping
   only SHA-256 and counts in memory.
2. Validate structure, count, **exact UID set**, and SHA-256 of each batch in a single
   pass. The exact-UID check is stronger than the current count-only check and is what
   makes F6 diagnosable rather than fatal.
3. Consolidate with `iterparse`, writing each `GBSeq` directly to the final temporary
   file. Never hold the whole tree.
4. Compute SHA-256 incrementally during the write.
5. Collapse passes 3, 4, 5 and 7, 8 into one streaming validator.
6. Return the already-validated snapshot after publication instead of reloading it.
7. Use the manifest's validated SHA-256 in notebook 01 instead of rehashing 363 MB.

Targets F4 and F5. Expected effect: peak memory from ~3 GB to well under 100 MB;
seven full passes reduced to one or two.

> Ordering note: P1 is placed before P3 because (a) it fixes a stability defect, not
> only a speed one, (b) larger batches make F4 worse by raising each intermediate
> parse from 0.88 MB to 4.37 MB of text, and (c) it costs no NCBI quota. If P0
> instrumentation shows local time is under 10% of the total, P1 and P3 may be
> swapped without consequence.

### P2. UID path

1. Run ESearch once with `usehistory=y`.
2. Retrieve UIDs from the frozen History set via EFetch with `rettype=uilist`,
   `retmode=text`, `retstart`, and `retmax=10000` (the documented ESearch/EFetch
   ceiling per NBK25499).
3. Requests: **43 -> 6**, a 86.05% reduction. Regulatory sleep: 5.46 s -> 0.65 s.
4. Apply the same timeout, deadline, and single retry layer already present on the
   XML path (**F11**). The UID path currently relies partly on Biopython's global
   `max_tries` / `sleep_between_tries`, which `_open_ncbi_entrez_efetch_once` exists
   specifically to bypass on the XML side.

Primary justification is **F1**, snapshot isolation. The time saving is secondary.

### P3. Larger XML batches

1. Move `batch_size` from 100 to 500. Keep execution sequential to isolate the effect.
2. Requests: **414 -> 83**, a 79.95% reduction. Regulatory sleep: 53.82 s -> 10.79 s.
3. `_open_ncbi_entrez_efetch_once` already switches to POST at `>= 200` UIDs, so no
   transport change is needed.
4. Hard ceiling: EFetch `retmax` saturates at 10,000. `batch_size` must never exceed
   it. Section 1.4 shows 500 is well inside safe payload territory.

**Deferred variant:** EPost of the frozen list followed by 83 History-backed EFetch
calls. `_validate_xml_record_uids` checks the exact UID set, but consolidated record
order follows batch arrival order. History set ordering is presumably submission
order, but this is **not documented**. If it diverges, `xml_file_sha256` changes with
no existing validation detecting it. Treat EPost as an experiment gated on an explicit
order check, not as an equal-weight candidate. Direct POST with 500 UIDs carries no
such risk.

### P4. Selective resume

1. Reintroduce the batch workspace as its own scope.
2. Write `batch_plan.json` containing the UID snapshot hash, batch size, ranges,
   request policy, and format version.
3. Reuse only batches whose hash, count, and exact UID set validate.
4. Selectively re-fetch missing, corrupt, or incompatible batches.
5. Purge incomplete artifacts before final publication.

Does not reduce the ideal-run time. Reduces expected time `E[T]` substantially given
the 5-retry and circuit-breaker regime, especially for a failure near batch 80 of 83.

### P5. Bounded concurrency

Blocked by: P0 item 4 (F3), P1 (ordered assembly), P4 (recoverability).

1. Add a global limiter keyed on request **start** instants, shared across workers.
   Do not use per-worker `time.sleep`.
2. Start at 2 concurrent transfers.
3. Move to 4 only if measurements show idle latency and zero 429s.
4. Conservative ceiling: 8 starts/s with an API key, 2.5/s without. NCBI's limit
   applies to start frequency and aggregates all activity under the same key, so
   concurrency within that budget is permitted.
5. Apply adaptive reduction on 429, 503, timeout, or circuit-breaker opening.
6. Test HTTP connection reuse with one session per worker (**F7**).

Concurrency-specific hazards to respect:

- `_active_ncbi_ssl_context` and `_active_ncbi_network_timeout_seconds` are
  `ContextVar`s. Executor threads **do not inherit** the caller's context.
  `copy_context()` must be taken at submit time, inside the
  `with _configured_ncbi_entrez_urlopen(...)` block. Omitting this silently discards
  a custom CA on TLS-intercepting networks.
- Results must be written into a list pre-sized and indexed by `batch_index`, never
  appended. Completion order must not become consolidation order, or `xml_file_sha256`
  stops being reproducible.
- `NCBIXmlCircuitBreaker` is already `RLock`-protected and shared. No change needed.
- `_fetch_ncbi_xml_attempt_with_deadline` spawns a watchdog thread per attempt, giving
  2P threads. Acceptable. The watchdog remains necessary even with a read timeout: a
  slow-but-alive stream defeats `SO_RCVTIMEO`.

---

## 5. Acceptance criteria

1. UIDs retrieved from the same History set are exactly equal to the frozen set.
2. Requested UIDs equal the UIDs extracted from the returned XML records, exactly.
3. Manifests, hashes, upstream provenance, immutable snapshots, and atomic
   publication are preserved.
4. Existing timeout, deadline, retry, partial-response, and circuit-breaker controls
   are preserved.
5. Zero 429 responses in the full trial.
6. Demonstrated reduction in total time, peak memory, and number of full passes over
   the XML.
7. `sha256(A) == sha256(B)` in the P0 `rettype` invariance test.
8. Explicit compatibility with prior manifests, or a documented format version bump.
9. `python -m unittest discover -s tests -q` passes, including
   `tests/test_ncbi_xml_snapshot.py`.

---

## 6. Measurement protocol

1. Use the existing immutable UID snapshot for all XML trials. Do not re-fetch UIDs
   to test XML changes.
2. Run a ~5,000 UID trial first, then a single full confirmation run.
3. Run the full confirmation inside NCBI's recommended window for large jobs:
   weekends, or 21:00-05:00 US Eastern on weekdays.
4. Hold one variable per trial. P3 in particular must be measured sequentially before
   P5 is enabled, or batch-size and concurrency effects become inseparable.

---

## 7. Risk register

| Risk | Likelihood | Impact | Mitigation |
| --- | --- | --- | --- |
| `rettype=gp` alters the payload | Low | High | P0 invariance test before any full run |
| Concurrency triggers F3 and kills a long run | High if F3 unfixed | High | F3 is a hard P0 gate for P5 |
| Unordered assembly changes `xml_file_sha256` | Medium | High | Index-addressed assembly; hash comparison in acceptance |
| ContextVar not propagated to workers | Medium | High on restricted networks | Explicit `copy_context()` at submit; test with a custom CA |
| EPost order differs from submission order | Unknown | High | Deferred variant; explicit order check required |
| Batch of 500 amplifies F6 blast radius | Medium | Medium | P1 exact-UID validation plus P4 resume, both before P3 |
| IP blocked for policy violation | Low | Severe | Conservative ceilings; registered `tool`/`email`; off-peak window |

---

## 8. Open decisions

1. **Manifest compatibility.** Changing `batch_size` changes `batch_size`,
   `batch_count`, and `batches[]` in the manifest. Consolidated XML content and order
   are unaffected if assembly is index-ordered. Decide whether to bump
   `snapshot_format_version` or to treat these fields as free.
2. **New HTTP dependency.** `httpx` or `requests` would provide keep-alive, separate
   connect/read timeouts, and native `verify=`, likely allowing removal of
   `_ensure_ncbi_entrez_urlopen_hook_installed` and both SSL/timeout `ContextVar`s.
   Net effect on line count is probably negative. Declining this forfeits F7 (~9%)
   but leaves P2, P3, and P5 intact.
3. **Scientific refresh.** Independent of performance: whether to re-run the query
   against current NCBI content. Per F12 this is not forced by the tree state.
4. **`GBSeq_references`.** 20.6% of payload with no identified consumer, and no API
   parameter to suppress it. Recorded for completeness; no action available.

---

## 9. Expected outcome

| Stage | Current | After P2/P3 | After P5 (P=4) |
| --- | --- | --- | --- |
| UID requests | 43 | 6 (-86.05%) | 6 |
| XML requests | 414 | 83 (-79.95%) | 83 |
| UID sleep | 5.46 s | 0.65 s | 0.65 s |
| XML sleep | 53.82 s | 10.79 s | 10.79 s |
| Consolidation peak memory | ~2.5-3.5 GB | <100 MB (P1) | <100 MB |
| Full XML passes | 7 + 1 copy | 1-2 (P1) | 1-2 |

Absolute wall-clock reduction cannot be stated until P0 supplies `alpha` and
`t_conn`. Based on the model in section 2, the remote stage should fall by roughly
20-25% from P3 and P7 alone, and by roughly 70-80% once bounded concurrency is
enabled. These are model projections, not measurements.

---

## 10. Recommended order

```
P0  instrumentation + F3 + F2 invariance test + F10       (2 NCBI requests)
P1  streaming persistence and validation                  (0 NCBI requests)
P2  UID History path + F11                                (6 NCBI requests)
P3  batch_size 100 -> 500, sequential                     (83 NCBI requests)
P4  selective resume                                      (0 NCBI requests)
P5  bounded concurrency + connection reuse                (83 NCBI requests)
```
