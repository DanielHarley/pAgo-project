import contextlib
import http.client
import io
import ssl
import sys
import tempfile
import threading
import time
import types
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError

if "Bio" not in sys.modules:
    fake_bio_module = types.ModuleType("Bio")
    fake_bio_module.__version__ = "test"
    fake_entrez_module = types.ModuleType("Bio.Entrez")
    fake_entrez_module.urlopen = lambda *args, **kwargs: None
    fake_seqio_module = types.ModuleType("Bio.SeqIO")
    fake_bio_module.Entrez = fake_entrez_module
    fake_bio_module.SeqIO = fake_seqio_module
    sys.modules["Bio"] = fake_bio_module
    sys.modules["Bio.Entrez"] = fake_entrez_module
    sys.modules["Bio.SeqIO"] = fake_seqio_module

import src.pago_pipeline.ncbi_api as ncbi_api
from src.pago_pipeline.ncbi_api import (
    NCBI_PROTEIN_GBSEQ_XML_RETTYPE,
    _build_xml_batch_context,
    _format_duration_seconds,
    NCBIRetrievalProgressReporter,
    _payload_indicates_ncbi_rate_limit,
    fetch_ncbi_protein_xml_batches,
    read_ncbi_xml_batch_payload_bytes,
)
from src.pago_pipeline.ncbi_telemetry import (
    NCBIAdaptiveConcurrencyGovernor,
    NCBIRequestStartRateLimiter,
)


class FakeFetchHandle:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.closed = False
        self.url = "https://ncbi.example.test/efetch"
        self.headers = {"content-type": "text/xml"}

    def read(self):
        return self.payload

    def close(self) -> None:
        self.closed = True


def build_gbset_payload(protein_uids: list[str]) -> bytes:
    gbseq_records = b"".join(
        f"<GBSeq><GBSeqid>gi|{protein_uid}</GBSeqid></GBSeq>".encode("utf-8")
        for protein_uid in protein_uids
    )
    return b"<GBSet>" + gbseq_records + b"</GBSet>"


class XmlFetchTestCase(unittest.TestCase):
    def setUp(self) -> None:
        with ncbi_api._default_ncbi_xml_circuit_breakers_lock:
            ncbi_api._default_ncbi_xml_circuit_breakers.clear()

        batch_spill_directory = tempfile.TemporaryDirectory()
        self.addCleanup(batch_spill_directory.cleanup)
        self.batch_spill_directory_name = batch_spill_directory.name

    def fetch_xml_batches(self, **fetch_parameters):
        fetch_parameters.setdefault("ncbi_email", "test@example.org")
        fetch_parameters.setdefault("ncbi_api_key", None)
        fetch_parameters.setdefault("request_delay_seconds", 0.0)
        fetch_parameters.setdefault(
            "batch_spill_parent_directory",
            self.batch_spill_directory_name,
        )
        return fetch_ncbi_protein_xml_batches(**fetch_parameters)


class EfetchRequestFormatTests(XmlFetchTestCase):
    def _capture_request_url(self, **fetch_parameters) -> str:
        captured_request_urls: list[str] = []

        def fake_urlopen(request, *args, **kwargs):
            captured_request_urls.append(request.full_url)
            return FakeFetchHandle(build_gbset_payload(["1001"]))

        with patch.object(ncbi_api.Entrez, "urlopen", fake_urlopen, create=True):
            with patch("src.pago_pipeline.ncbi_api.random.uniform", return_value=0.0):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    self.fetch_xml_batches(
                        protein_uids=["1001"],
                        batch_size=1,
                        **fetch_parameters,
                    )

        return captured_request_urls[0]

    def test_gbseq_xml_rettype_is_requested_explicitly(self) -> None:
        request_url = self._capture_request_url()

        self.assertIn(f"rettype={NCBI_PROTEIN_GBSEQ_XML_RETTYPE}", request_url)
        self.assertIn("retmode=xml", request_url)
        self.assertIn("db=protein", request_url)

    def test_rettype_can_be_omitted_to_reproduce_the_previous_default(self) -> None:
        request_url = self._capture_request_url(rettype=None)

        self.assertNotIn("rettype=", request_url)
        self.assertIn("retmode=xml", request_url)

    def test_manifest_facing_result_records_the_requested_rettype(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            return_value=FakeFetchHandle(build_gbset_payload(["1001"])),
        ):
            with patch("src.pago_pipeline.ncbi_api.random.uniform", return_value=0.0):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_result = self.fetch_xml_batches(
                        protein_uids=["1001"],
                        batch_size=1,
                    )

        self.assertEqual(fetch_result.rettype, NCBI_PROTEIN_GBSEQ_XML_RETTYPE)
        self.assertEqual(fetch_result.retmode, "xml")


class RateLimitResponseTests(XmlFetchTestCase):
    def test_rate_limit_body_is_detected(self) -> None:
        self.assertTrue(
            _payload_indicates_ncbi_rate_limit(
                payload_bytes=b'{"error":"API rate limit exceeded","count":"11"}',
            )
        )
        self.assertTrue(
            _payload_indicates_ncbi_rate_limit(
                payload_bytes=b'{"error":"Request rate limit reached"}',
            )
        )
        self.assertFalse(
            _payload_indicates_ncbi_rate_limit(
                payload_bytes=build_gbset_payload(["1001"]),
            )
        )

    def test_rate_limit_body_is_retried_with_an_extended_backoff(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                FakeFetchHandle(b'{"error":"API rate limit exceeded","count":"11"}'),
                FakeFetchHandle(build_gbset_payload(["1001"])),
            ],
        ) as efetch_mock:
            with patch("src.pago_pipeline.ncbi_api.random.uniform", return_value=0.0):
                with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                    fetch_result = self.fetch_xml_batches(
                        protein_uids=["1001"],
                        batch_size=1,
                        max_retry_attempts=2,
                        retry_backoff_initial_seconds=0.25,
                        rate_limit_backoff_seconds=5.0,
                        batch_deadline_seconds=60.0,
                    )

        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(fetch_result.batch_count, 1)
        sleep_mock.assert_any_call(5.0)

    def test_rate_limit_failures_are_counted_in_telemetry(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                FakeFetchHandle(b'{"error":"API rate limit exceeded","count":"11"}'),
                FakeFetchHandle(build_gbset_payload(["1001"])),
            ],
        ):
            with patch("src.pago_pipeline.ncbi_api.random.uniform", return_value=0.0):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_result = self.fetch_xml_batches(
                        protein_uids=["1001"],
                        batch_size=1,
                        max_retry_attempts=2,
                        retry_backoff_initial_seconds=0.25,
                        batch_deadline_seconds=60.0,
                    )

        xml_stage_telemetry = fetch_result.retrieval_telemetry["stages"][
            "xml_batches"
        ]
        self.assertEqual(xml_stage_telemetry["failure_counts"]["rate_limit"], 1)
        self.assertEqual(xml_stage_telemetry["retry_count"], 1)
        self.assertEqual(xml_stage_telemetry["request_count"], 2)


class TruncatedResponseTests(XmlFetchTestCase):
    """
    A response that ends early must be retried, not treated as permanent.

    The request was accepted and the server began answering it, so nothing
    about it is unrepeatable. Larger batches make this far more likely: a
    multi-megabyte body spends much longer exposed to a mid-stream drop than a
    small one.
    """

    def _assert_error_is_retried(self, transport_error: Exception) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                transport_error,
                FakeFetchHandle(build_gbset_payload(["1001"])),
            ],
        ) as efetch_mock:
            with patch("src.pago_pipeline.ncbi_api.random.uniform", return_value=0.0):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_result = self.fetch_xml_batches(
                        protein_uids=["1001"],
                        batch_size=1,
                        max_retry_attempts=2,
                        retry_backoff_initial_seconds=0.1,
                    )

        self.assertEqual(efetch_mock.call_count, 2)
        self.assertEqual(fetch_result.batch_count, 1)
        return fetch_result

    def test_incomplete_read_is_retried(self) -> None:
        self._assert_error_is_retried(
            http.client.IncompleteRead(b"x" * 4_112_725, 700_000)
        )

    def test_remote_disconnected_is_retried(self) -> None:
        self._assert_error_is_retried(
            http.client.RemoteDisconnected(
                "Remote end closed connection without response"
            )
        )

    def test_connection_reset_is_retried(self) -> None:
        self._assert_error_is_retried(ConnectionResetError(10054, "forcibly closed"))

    def test_mid_stream_tls_eof_is_retried_not_reported_as_a_ca_problem(
        self,
    ) -> None:
        self._assert_error_is_retried(
            ssl.SSLEOFError("EOF occurred in violation of protocol")
        )

    def test_certificate_verification_failure_is_still_a_configuration_error(
        self,
    ) -> None:
        certificate_error = ssl.SSLCertVerificationError(
            "certificate verify failed: unable to get local issuer certificate"
        )

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=certificate_error,
        ) as efetch_mock:
            with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                with self.assertRaisesRegex(
                    RuntimeError,
                    "NCBI_SSL_CERT_FILE / NCBI_SSL_CERT_DIR",
                ):
                    self.fetch_xml_batches(
                        protein_uids=["1001"],
                        batch_size=1,
                        max_retry_attempts=3,
                    )

        efetch_mock.assert_called_once()

    def test_truncated_responses_are_counted_separately_in_telemetry(self) -> None:
        fetch_result = self._assert_error_is_retried(
            http.client.IncompleteRead(b"x" * 100, 50)
        )

        xml_stage_telemetry = fetch_result.retrieval_telemetry["stages"][
            "xml_batches"
        ]
        self.assertEqual(
            xml_stage_telemetry["failure_counts"]["truncated_response"],
            1,
        )
        self.assertEqual(xml_stage_telemetry["failure_counts"]["other"], 0)


class BatchContextTests(unittest.TestCase):
    def test_batch_context_summarizes_the_uid_range_by_default(self) -> None:
        batch_context = _build_xml_batch_context(
            batch_index=2,
            total_batch_count=10,
            batch_start_index=500,
            current_batch_protein_uids=[str(uid) for uid in range(1000, 1500)],
        )

        self.assertIn("protein_uids=1000..1499 (n=500)", batch_context)
        self.assertLess(len(batch_context), 200)

    def test_batch_context_can_include_the_full_uid_list(self) -> None:
        batch_context = _build_xml_batch_context(
            batch_index=1,
            total_batch_count=1,
            batch_start_index=0,
            current_batch_protein_uids=["1001", "1002"],
            include_full_protein_uid_list=True,
        )

        self.assertIn("protein_uids=['1001', '1002']", batch_context)


class BatchProgressReportingTests(XmlFetchTestCase):
    def _build_reporter(self, *, total_batch_count: int):
        self.fake_clock = {"value": 0.0}
        return NCBIRetrievalProgressReporter(
            total_unit_count=total_batch_count,
            inline_updates_enabled=False,
            monotonic_function=lambda: self.fake_clock["value"],
        )

    def _capture_progress_line(self, progress_reporter, **report_parameters) -> str:
        captured_output = io.StringIO()
        with contextlib.redirect_stdout(captured_output):
            progress_reporter.report_completed_units(**report_parameters)

        return captured_output.getvalue().strip()

    def test_durations_are_formatted_for_a_reader(self) -> None:
        self.assertEqual(_format_duration_seconds(4.567), "4.6s")
        self.assertEqual(_format_duration_seconds(59.9), "59.9s")
        self.assertEqual(_format_duration_seconds(75), "1m15s")
        self.assertEqual(_format_duration_seconds(3675), "1h01m15s")

    def test_progress_line_reports_completion_throughput_and_finish_time(
        self,
    ) -> None:
        progress_reporter = self._build_reporter(total_batch_count=4)

        self.fake_clock["value"] = 8.0
        progress_line = self._capture_progress_line(
            progress_reporter,
            completed_unit_delta=1,
            item_label="#1",
            response_byte_count=2 * 1024 * 1024,
            round_trip_latency_seconds=3.25,
        )

        self.assertIn("1/4  25.0%", progress_line)
        self.assertIn("#1 3.2s 2.00 MB", progress_line)
        self.assertIn("elapsed", progress_line)
        self.assertIn("2.0 MB at", progress_line)
        self.assertIn("left", progress_line)
        self.assertIn("ends ", progress_line)

    def test_reused_batches_do_not_drive_the_projection(self) -> None:
        progress_reporter = self._build_reporter(total_batch_count=4)

        self.fake_clock["value"] = 4.0
        reused_batch_line = self._capture_progress_line(
            progress_reporter,
            completed_unit_delta=1,
            item_label="#1",
            response_byte_count=4096,
            required_a_request=False,
        )

        # A resumed run completes reused batches instantly. Counting them in the
        # rate would project a finish time that no real request can meet.
        self.assertIn("1/4  25.0%", reused_batch_line)
        self.assertIn("#1 reused", reused_batch_line)
        self.assertNotIn("ends ", reused_batch_line)
        self.assertEqual(progress_reporter.network_byte_count, 0)

        self.fake_clock["value"] = 8.0
        fetched_batch_line = self._capture_progress_line(
            progress_reporter,
            completed_unit_delta=1,
            item_label="#2",
            response_byte_count=4096,
        )

        self.assertIn("2/4  50.0%", fetched_batch_line)
        self.assertIn("ends ", fetched_batch_line)
        self.assertEqual(progress_reporter.network_byte_count, 4096)

    def test_last_batch_omits_the_projection(self) -> None:
        progress_reporter = self._build_reporter(total_batch_count=1)

        self.fake_clock["value"] = 5.0
        progress_line = self._capture_progress_line(
            progress_reporter,
            completed_unit_delta=1,
            item_label="#1",
            response_byte_count=1024,
        )

        self.assertIn("1/1 100.0%", progress_line)
        self.assertNotIn("left", progress_line)
        self.assertNotIn("ends ", progress_line)

    def test_fetch_run_emits_progress_and_a_total_elapsed_summary(self) -> None:
        captured_output = io.StringIO()

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                FakeFetchHandle(build_gbset_payload(["1001", "1002"])),
                FakeFetchHandle(build_gbset_payload(["1003", "1004"])),
            ],
        ):
            with patch("src.pago_pipeline.ncbi_api.random.uniform", return_value=0.0):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    with contextlib.redirect_stdout(captured_output):
                        self.fetch_xml_batches(
                            protein_uids=["1001", "1002", "1003", "1004"],
                            batch_size=2,
                        )

        retrieval_log = captured_output.getvalue()
        self.assertIn("1/2  50.0%", retrieval_log)
        self.assertIn("2/2 100.0%", retrieval_log)
        # The per-request chatter is behind the verbose flag now.
        self.assertNotIn("Starting XML request", retrieval_log)
        self.assertNotIn("Received XML response", retrieval_log)
        self.assertRegex(
            retrieval_log,
            r"Fetched 2 XML batches for 4 protein UIDs in [\d.]+s\.",
        )
        self.assertIn("Transferred", retrieval_log)


class BatchSizeCeilingTests(XmlFetchTestCase):
    def test_batch_size_above_the_efetch_ceiling_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "EFetch paging ceiling"):
            self.fetch_xml_batches(
                protein_uids=["1001"],
                batch_size=10_001,
            )

    def test_default_batch_size_is_five_hundred(self) -> None:
        import inspect

        signature = inspect.signature(fetch_ncbi_protein_xml_batches)
        self.assertEqual(signature.parameters["batch_size"].default, 500)


class BoundedConcurrencyTests(XmlFetchTestCase):
    def test_results_follow_plan_order_when_batches_complete_out_of_order(
        self,
    ) -> None:
        protein_uids = [str(1000 + index) for index in range(6)]
        completion_order: list[int] = []
        completion_lock = threading.Lock()

        def fake_efetch(**parameters):
            requested_protein_uids = list(parameters["protein_uids"])
            # The first batch is deliberately the slowest so that completion
            # order cannot coincide with plan order.
            if requested_protein_uids[0] == protein_uids[0]:
                time.sleep(0.15)

            with completion_lock:
                completion_order.append(int(requested_protein_uids[0]))

            return FakeFetchHandle(build_gbset_payload(requested_protein_uids))

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            fake_efetch,
        ):
            fetch_result = self.fetch_xml_batches(
                protein_uids=protein_uids,
                batch_size=2,
                max_concurrent_requests=3,
                max_request_starts_per_second=1000.0,
            )

        self.assertEqual(
            [xml_batch.batch_index for xml_batch in fetch_result.xml_batches],
            [1, 2, 3],
        )
        self.assertEqual(
            [xml_batch.protein_uids for xml_batch in fetch_result.xml_batches],
            [["1000", "1001"], ["1002", "1003"], ["1004", "1005"]],
        )
        self.assertNotEqual(completion_order[0], 1000)
        self.assertEqual(
            read_ncbi_xml_batch_payload_bytes(
                xml_batch=fetch_result.xml_batches[0],
            ),
            build_gbset_payload(["1000", "1001"]),
        )

    def test_worker_threads_inherit_the_configured_ssl_context(self) -> None:
        observed_ssl_contexts: list[object] = []
        ssl_context_sentinel = object()

        def fake_efetch(**parameters):
            observed_ssl_contexts.append(
                ncbi_api._active_ncbi_ssl_context.get()
            )
            return FakeFetchHandle(
                build_gbset_payload(list(parameters["protein_uids"]))
            )

        with tempfile.TemporaryDirectory() as temporary_directory_name:
            ssl_ca_file_path = Path(temporary_directory_name) / "corporate-ca.pem"
            ssl_ca_file_path.write_text("dummy-ca", encoding="utf-8")

            with patch(
                "src.pago_pipeline.ncbi_api.ssl.create_default_context",
                return_value=ssl_context_sentinel,
            ):
                with patch(
                    "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
                    fake_efetch,
                ):
                    self.fetch_xml_batches(
                        protein_uids=["1001", "1002", "1003", "1004"],
                        batch_size=1,
                        max_concurrent_requests=2,
                        max_request_starts_per_second=1000.0,
                        ssl_ca_file=str(ssl_ca_file_path),
                    )

        self.assertEqual(len(observed_ssl_contexts), 4)
        for observed_ssl_context in observed_ssl_contexts:
            self.assertIs(observed_ssl_context, ssl_context_sentinel)

    def test_a_failing_batch_aborts_the_concurrent_run(self) -> None:
        permanent_http_error = HTTPError(
            "https://ncbi.example.test/efetch",
            400,
            "Bad Request",
            None,
            None,
        )
        self.addCleanup(permanent_http_error.close)

        def fake_efetch(**parameters):
            requested_protein_uids = list(parameters["protein_uids"])
            if requested_protein_uids == ["1003"]:
                raise permanent_http_error
            return FakeFetchHandle(build_gbset_payload(requested_protein_uids))

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            fake_efetch,
        ):
            with self.assertRaisesRegex(
                RuntimeError,
                "Permanent NCBI XML batch failure",
            ):
                self.fetch_xml_batches(
                    protein_uids=["1001", "1002", "1003", "1004"],
                    batch_size=1,
                    max_concurrent_requests=2,
                    max_request_starts_per_second=1000.0,
                )

    def test_concurrent_runs_replace_the_per_request_sleep_with_a_limiter(
        self,
    ) -> None:
        def fake_efetch(**parameters):
            return FakeFetchHandle(
                build_gbset_payload(list(parameters["protein_uids"]))
            )

        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            fake_efetch,
        ):
            with patch("src.pago_pipeline.ncbi_api.time.sleep") as sleep_mock:
                fetch_result = self.fetch_xml_batches(
                    protein_uids=["1001", "1002"],
                    batch_size=1,
                    max_concurrent_requests=2,
                    max_request_starts_per_second=1000.0,
                )

        sleep_mock.assert_not_called()
        self.assertEqual(fetch_result.max_concurrent_requests, 2)
        self.assertEqual(fetch_result.max_request_starts_per_second, 1000.0)


class RequestStartRateLimiterTests(unittest.TestCase):
    def test_request_starts_are_spaced_by_the_configured_interval(self) -> None:
        fake_clock = {"value": 0.0}
        observed_sleeps: list[float] = []

        def fake_sleep(sleep_seconds: float) -> None:
            observed_sleeps.append(sleep_seconds)
            fake_clock["value"] += sleep_seconds

        rate_limiter = NCBIRequestStartRateLimiter(
            max_request_starts_per_second=2.0,
            sleep_function=fake_sleep,
            monotonic_function=lambda: fake_clock["value"],
        )

        for _ in range(3):
            rate_limiter.acquire_request_start_slot()

        self.assertEqual(observed_sleeps, [0.5, 0.5])
        self.assertEqual(rate_limiter.acquired_start_count, 3)

    def test_rate_can_be_reduced_and_restored(self) -> None:
        rate_limiter = NCBIRequestStartRateLimiter(
            max_request_starts_per_second=8.0,
        )

        rate_limiter.set_max_request_starts_per_second(
            max_request_starts_per_second=2.0,
        )

        self.assertEqual(rate_limiter.max_request_starts_per_second, 2.0)


class AdaptiveConcurrencyGovernorTests(unittest.TestCase):
    def test_overload_reduces_concurrency_and_request_rate(self) -> None:
        rate_limiter = NCBIRequestStartRateLimiter(
            max_request_starts_per_second=4.0,
        )
        governor = NCBIAdaptiveConcurrencyGovernor(
            max_concurrent_requests=4,
            rate_limiter=rate_limiter,
            successes_before_recovery=2,
        )

        governor.record_overload_signal()

        self.assertEqual(governor.permitted_concurrent_requests, 3)
        self.assertEqual(rate_limiter.max_request_starts_per_second, 2.0)
        self.assertEqual(governor.reduction_event_count, 1)

    def test_recovery_requires_a_run_of_clean_responses(self) -> None:
        rate_limiter = NCBIRequestStartRateLimiter(
            max_request_starts_per_second=4.0,
        )
        governor = NCBIAdaptiveConcurrencyGovernor(
            max_concurrent_requests=4,
            rate_limiter=rate_limiter,
            successes_before_recovery=2,
        )
        governor.record_overload_signal()

        governor.record_successful_request()
        self.assertEqual(governor.permitted_concurrent_requests, 3)

        governor.record_successful_request()
        self.assertEqual(governor.permitted_concurrent_requests, 4)
        self.assertEqual(rate_limiter.max_request_starts_per_second, 4.0)

    def test_concurrency_never_drops_below_one(self) -> None:
        governor = NCBIAdaptiveConcurrencyGovernor(max_concurrent_requests=2)

        for _ in range(5):
            governor.record_overload_signal()

        self.assertEqual(governor.permitted_concurrent_requests, 1)

    def test_request_slots_bound_the_number_of_in_flight_requests(self) -> None:
        governor = NCBIAdaptiveConcurrencyGovernor(max_concurrent_requests=2)
        in_flight_counter = {"current": 0, "max": 0}
        counter_lock = threading.Lock()
        release_workers = threading.Event()

        def occupy_request_slot() -> None:
            with governor.request_slot():
                with counter_lock:
                    in_flight_counter["current"] += 1
                    in_flight_counter["max"] = max(
                        in_flight_counter["max"],
                        in_flight_counter["current"],
                    )
                release_workers.wait(timeout=1.0)
                with counter_lock:
                    in_flight_counter["current"] -= 1

        worker_threads = [
            threading.Thread(target=occupy_request_slot) for _ in range(4)
        ]
        for worker_thread in worker_threads:
            worker_thread.start()

        time.sleep(0.1)
        release_workers.set()
        for worker_thread in worker_threads:
            worker_thread.join(timeout=2.0)

        self.assertEqual(in_flight_counter["max"], 2)


class FakeHttpResponse:
    def __init__(self, *, status: int, payload: bytes) -> None:
        self.status = status
        self.reason = "OK" if status < 400 else "Error"
        self.headers = {"content-type": "text/xml"}
        self.payload = payload
        self.closed = False

    def read(self, *args, **kwargs):
        return self.payload

    def close(self) -> None:
        self.closed = True


class FakeHttpsConnection:
    created_connection_count = 0

    def __init__(self, host, port=None, timeout=None, context=None) -> None:
        FakeHttpsConnection.created_connection_count += 1
        self.host = host
        self.timeout = timeout
        self.context = context
        self.requests: list[tuple[str, str]] = []
        self.closed = False
        self.response_status = 200

    def connect(self) -> None:
        return None

    def request(self, method, target, body=None, headers=None) -> None:
        self.requests.append((method, target))

    def getresponse(self):
        return FakeHttpResponse(
            status=self.response_status,
            payload=b"<GBSet><GBSeq /></GBSet>",
        )

    def close(self) -> None:
        self.closed = True


class PersistentHttpsTransportTests(unittest.TestCase):
    def setUp(self) -> None:
        FakeHttpsConnection.created_connection_count = 0
        self.transport = ncbi_api._NCBIPersistentHttpsTransport()
        self.addCleanup(self.transport.discard_all_connections)

    def _open_request(self):
        return self.transport.open_request(
            request_url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein",
            http_method="GET",
            request_body=None,
            timeout_seconds=30.0,
            ssl_context=None,
        )

    def test_connection_is_reused_across_consecutive_requests(self) -> None:
        with patch.object(
            ncbi_api.http.client,
            "HTTPSConnection",
            FakeHttpsConnection,
        ):
            for _ in range(3):
                response = self._open_request()
                response.read()
                response.close()

        self.assertEqual(FakeHttpsConnection.created_connection_count, 1)

    def test_first_response_reports_the_connection_setup_cost(self) -> None:
        with patch.object(
            ncbi_api.http.client,
            "HTTPSConnection",
            FakeHttpsConnection,
        ):
            first_response = self._open_request()
            first_response.read()
            first_response.close()

            second_response = self._open_request()
            second_response.read()
            second_response.close()

        self.assertIsNotNone(first_response.ncbi_connect_seconds)
        self.assertIsNone(second_response.ncbi_connect_seconds)

    def test_all_connections_are_closed_including_other_threads(self) -> None:
        worker_connection_opened = threading.Event()

        def open_request_on_worker_thread() -> None:
            response = self._open_request()
            response.read()
            response.close()
            worker_connection_opened.set()

        with patch.object(
            ncbi_api.http.client,
            "HTTPSConnection",
            FakeHttpsConnection,
        ):
            worker_thread = threading.Thread(target=open_request_on_worker_thread)
            worker_thread.start()
            worker_thread.join(timeout=2.0)

            owning_thread_response = self._open_request()
            owning_thread_response.read()
            owning_thread_response.close()

            opened_connections = list(self.transport._open_connections)
            self.transport.discard_all_connections()

        self.assertTrue(worker_connection_opened.is_set())
        self.assertEqual(len(opened_connections), 2)
        for opened_connection in opened_connections:
            self.assertTrue(opened_connection.closed)

    def test_abandoned_response_body_discards_the_connection(self) -> None:
        with patch.object(
            ncbi_api.http.client,
            "HTTPSConnection",
            FakeHttpsConnection,
        ):
            first_response = self._open_request()
            first_response.close()

            second_response = self._open_request()
            second_response.read()
            second_response.close()

        self.assertEqual(FakeHttpsConnection.created_connection_count, 2)

    def test_error_status_is_raised_as_an_http_error(self) -> None:
        class FailingConnection(FakeHttpsConnection):
            def getresponse(self):
                return FakeHttpResponse(status=429, payload=b"")

        with patch.object(
            ncbi_api.http.client,
            "HTTPSConnection",
            FailingConnection,
        ):
            with self.assertRaises(HTTPError) as raised_error:
                self._open_request()

        self.addCleanup(raised_error.exception.close)
        self.assertEqual(raised_error.exception.code, 429)

    def test_efetch_can_route_through_the_persistent_transport(self) -> None:
        with patch.object(
            ncbi_api._ncbi_persistent_https_transport,
            "open_request",
        ) as open_request_mock:
            ncbi_api._open_ncbi_entrez_efetch_once(
                database_name="protein",
                retmode="xml",
                rettype="gp",
                protein_uids=["1001"],
                reuse_http_connection=True,
            )

        open_request_mock.assert_called_once()
        self.assertEqual(
            open_request_mock.call_args.kwargs["http_method"],
            "GET",
        )
        self.assertIn(
            "rettype=gp",
            open_request_mock.call_args.kwargs["request_url"],
        )


class RetrievalTelemetryTests(XmlFetchTestCase):
    def test_successful_run_reports_request_volume_and_latency(self) -> None:
        with patch(
            "src.pago_pipeline.ncbi_api._open_ncbi_entrez_efetch_once",
            side_effect=[
                FakeFetchHandle(build_gbset_payload(["1001", "1002"])),
                FakeFetchHandle(build_gbset_payload(["1003", "1004"])),
            ],
        ):
            with patch("src.pago_pipeline.ncbi_api.random.uniform", return_value=0.0):
                with patch("src.pago_pipeline.ncbi_api.time.sleep"):
                    fetch_result = self.fetch_xml_batches(
                        protein_uids=["1001", "1002", "1003", "1004"],
                        batch_size=2,
                        request_delay_seconds=0.1,
                    )

        xml_stage_telemetry = fetch_result.retrieval_telemetry["stages"][
            "xml_batches"
        ]
        self.assertEqual(xml_stage_telemetry["request_count"], 2)
        self.assertEqual(xml_stage_telemetry["retry_count"], 0)
        self.assertGreater(xml_stage_telemetry["response_byte_total"], 0)
        self.assertAlmostEqual(
            xml_stage_telemetry["sleep_seconds_total"],
            0.2,
            places=6,
        )
        self.assertEqual(
            xml_stage_telemetry["time_to_first_byte_latency"]["sample_count"],
            2,
        )
        self.assertEqual(
            xml_stage_telemetry["response_read_latency"]["sample_count"],
            2,
        )
        self.assertIn("process", fetch_result.retrieval_telemetry)


if __name__ == "__main__":
    unittest.main()
