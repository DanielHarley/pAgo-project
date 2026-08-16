from __future__ import annotations

import statistics
import threading
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

try:  # pragma: no cover - psutil is an optional runtime dependency
    import psutil as _psutil
except Exception:  # pragma: no cover - telemetry must never break retrieval
    _psutil = None


TELEMETRY_FORMAT_VERSION = "1.0"

UID_SEARCH_STAGE_NAME = "uid_initial_search"
UID_PAGE_STAGE_NAME = "uid_pages"
XML_BATCH_STAGE_NAME = "xml_batches"

RATE_LIMIT_FAILURE_KIND = "rate_limit"
HTTP_429_FAILURE_KIND = "http_429"
HTTP_5XX_FAILURE_KIND = "http_5xx"
TIMEOUT_FAILURE_KIND = "timeout"
TRUNCATED_RESPONSE_FAILURE_KIND = "truncated_response"
DEADLINE_FAILURE_KIND = "deadline"
VALIDATION_FAILURE_KIND = "response_validation"
OTHER_FAILURE_KIND = "other"

_FAILURE_KINDS = (
    RATE_LIMIT_FAILURE_KIND,
    HTTP_429_FAILURE_KIND,
    HTTP_5XX_FAILURE_KIND,
    TIMEOUT_FAILURE_KIND,
    TRUNCATED_RESPONSE_FAILURE_KIND,
    DEADLINE_FAILURE_KIND,
    VALIDATION_FAILURE_KIND,
    OTHER_FAILURE_KIND,
)


def _percentile_of_sorted_samples(
    *,
    sorted_samples: List[float],
    percentile_fraction: float,
) -> Optional[float]:
    """
    Return a nearest-rank percentile from an already sorted sample list.
    """
    if not sorted_samples:
        return None

    sample_count = len(sorted_samples)
    nearest_rank = int(percentile_fraction * sample_count + 0.999999)
    nearest_rank = max(1, min(sample_count, nearest_rank))
    return sorted_samples[nearest_rank - 1]


def _summarize_duration_samples(
    *,
    duration_samples: List[float],
) -> Dict[str, Any]:
    """
    Summarize a duration sample list without pulling in a numeric dependency.
    """
    if not duration_samples:
        return {
            "sample_count": 0,
            "total_seconds": 0.0,
            "median_seconds": None,
            "p95_seconds": None,
            "max_seconds": None,
        }

    sorted_samples = sorted(duration_samples)
    return {
        "sample_count": len(sorted_samples),
        "total_seconds": round(sum(sorted_samples), 6),
        "median_seconds": round(statistics.median(sorted_samples), 6),
        "p95_seconds": round(
            _percentile_of_sorted_samples(
                sorted_samples=sorted_samples,
                percentile_fraction=0.95,
            ),
            6,
        ),
        "max_seconds": round(sorted_samples[-1], 6),
    }


@dataclass
class _StageAccumulator:
    """
    Mutable per-stage counters accumulated while a retrieval stage runs.
    """
    stage_name: str
    request_count: int = 0
    retry_count: int = 0
    reused_batch_count: int = 0
    response_byte_total: int = 0
    sleep_seconds_total: float = 0.0
    wall_seconds_total: float = 0.0
    total_latency_samples: List[float] = field(default_factory=list)
    connect_latency_samples: List[float] = field(default_factory=list)
    ttfb_latency_samples: List[float] = field(default_factory=list)
    read_latency_samples: List[float] = field(default_factory=list)
    failure_counts: Dict[str, int] = field(
        default_factory=lambda: {kind: 0 for kind in _FAILURE_KINDS}
    )

    def build_summary(self) -> Dict[str, Any]:
        return {
            "request_count": self.request_count,
            "retry_count": self.retry_count,
            "reused_batch_count": self.reused_batch_count,
            "response_byte_total": self.response_byte_total,
            "sleep_seconds_total": round(self.sleep_seconds_total, 6),
            "wall_seconds_total": round(self.wall_seconds_total, 6),
            "total_latency": _summarize_duration_samples(
                duration_samples=self.total_latency_samples,
            ),
            "connect_latency": _summarize_duration_samples(
                duration_samples=self.connect_latency_samples,
            ),
            "time_to_first_byte_latency": _summarize_duration_samples(
                duration_samples=self.ttfb_latency_samples,
            ),
            "response_read_latency": _summarize_duration_samples(
                duration_samples=self.read_latency_samples,
            ),
            "failure_counts": dict(self.failure_counts),
        }


class NCBIRetrievalTelemetry:
    """
    Thread-safe collector for NCBI retrieval timing, volume, and failure counters.

    The collector is deliberately passive: it never raises for unknown stages or
    unavailable operating-system counters, because instrumentation must not be
    able to fail a retrieval run it is only observing.
    """

    def __init__(self) -> None:
        self._state_lock = threading.RLock()
        self._stages: Dict[str, _StageAccumulator] = {}
        self._local_phase_seconds: Dict[str, float] = {}
        self._observed_max_rss_bytes: Optional[int] = None
        self._initial_disk_read_bytes: Optional[int] = None
        self._initial_disk_write_bytes: Optional[int] = None
        self._final_disk_read_bytes: Optional[int] = None
        self._final_disk_write_bytes: Optional[int] = None
        self._process_handle = self._resolve_process_handle()
        self._capture_initial_disk_counters()
        self.sample_resource_usage()

    @staticmethod
    def _resolve_process_handle():
        if _psutil is None:
            return None

        try:
            return _psutil.Process()
        except Exception:
            return None

    def _capture_initial_disk_counters(self) -> None:
        disk_counters = self._read_disk_counters()
        if disk_counters is None:
            return

        self._initial_disk_read_bytes, self._initial_disk_write_bytes = disk_counters

    def _read_disk_counters(self) -> Optional[tuple[int, int]]:
        if self._process_handle is None:
            return None

        try:
            io_counters = self._process_handle.io_counters()
        except Exception:
            return None

        read_bytes = getattr(io_counters, "read_bytes", None)
        write_bytes = getattr(io_counters, "write_bytes", None)
        if read_bytes is None or write_bytes is None:
            return None

        return int(read_bytes), int(write_bytes)

    def _stage(self, stage_name: str) -> _StageAccumulator:
        stage = self._stages.get(stage_name)
        if stage is None:
            stage = _StageAccumulator(stage_name=stage_name)
            self._stages[stage_name] = stage

        return stage

    def record_request(
        self,
        *,
        stage_name: str,
        total_seconds: float,
        response_byte_count: int = 0,
        connect_seconds: Optional[float] = None,
        time_to_first_byte_seconds: Optional[float] = None,
        response_read_seconds: Optional[float] = None,
    ) -> None:
        """
        Record one completed HTTP request attempt against NCBI.
        """
        with self._state_lock:
            stage = self._stage(stage_name)
            stage.request_count += 1
            stage.response_byte_total += int(response_byte_count)
            stage.total_latency_samples.append(float(total_seconds))

            if connect_seconds is not None:
                stage.connect_latency_samples.append(float(connect_seconds))
            if time_to_first_byte_seconds is not None:
                stage.ttfb_latency_samples.append(
                    float(time_to_first_byte_seconds)
                )
            if response_read_seconds is not None:
                stage.read_latency_samples.append(float(response_read_seconds))

    def record_retry(self, *, stage_name: str) -> None:
        with self._state_lock:
            self._stage(stage_name).retry_count += 1

    def record_reused_batch(self, *, stage_name: str) -> None:
        with self._state_lock:
            self._stage(stage_name).reused_batch_count += 1

    def record_failure(self, *, stage_name: str, failure_kind: str) -> None:
        resolved_failure_kind = (
            failure_kind if failure_kind in _FAILURE_KINDS else OTHER_FAILURE_KIND
        )
        with self._state_lock:
            stage = self._stage(stage_name)
            stage.failure_counts[resolved_failure_kind] += 1

    def record_sleep(self, *, stage_name: str, sleep_seconds: float) -> None:
        with self._state_lock:
            self._stage(stage_name).sleep_seconds_total += float(sleep_seconds)

    def record_stage_wall_seconds(
        self,
        *,
        stage_name: str,
        wall_seconds: float,
    ) -> None:
        with self._state_lock:
            self._stage(stage_name).wall_seconds_total += float(wall_seconds)

    def record_local_phase_seconds(
        self,
        *,
        phase_name: str,
        elapsed_seconds: float,
    ) -> None:
        with self._state_lock:
            self._local_phase_seconds[phase_name] = round(
                self._local_phase_seconds.get(phase_name, 0.0)
                + float(elapsed_seconds),
                6,
            )

    @contextmanager
    def measure_local_phase(self, *, phase_name: str):
        """
        Measure one local (non-network) processing phase such as consolidation.
        """
        started_at_seconds = time.perf_counter()
        try:
            yield
        finally:
            self.record_local_phase_seconds(
                phase_name=phase_name,
                elapsed_seconds=time.perf_counter() - started_at_seconds,
            )
            self.sample_resource_usage()

    @contextmanager
    def measure_stage(self, *, stage_name: str):
        """
        Measure the wall-clock duration of one complete retrieval stage.
        """
        started_at_seconds = time.perf_counter()
        try:
            yield
        finally:
            self.record_stage_wall_seconds(
                stage_name=stage_name,
                wall_seconds=time.perf_counter() - started_at_seconds,
            )
            self.sample_resource_usage()

    def sample_resource_usage(self) -> Optional[int]:
        """
        Sample resident memory at a checkpoint and retain the observed maximum.
        """
        if self._process_handle is None:
            return None

        try:
            resident_set_size_bytes = int(
                self._process_handle.memory_info().rss
            )
        except Exception:
            return None

        with self._state_lock:
            if (
                self._observed_max_rss_bytes is None
                or resident_set_size_bytes > self._observed_max_rss_bytes
            ):
                self._observed_max_rss_bytes = resident_set_size_bytes

        return resident_set_size_bytes

    def finalize(self) -> None:
        """
        Capture terminal process counters before a summary is built.
        """
        self.sample_resource_usage()
        disk_counters = self._read_disk_counters()
        if disk_counters is None:
            return

        with self._state_lock:
            self._final_disk_read_bytes, self._final_disk_write_bytes = disk_counters

    def _build_process_summary(self) -> Dict[str, Any]:
        disk_read_bytes: Optional[int] = None
        disk_write_bytes: Optional[int] = None

        if (
            self._initial_disk_read_bytes is not None
            and self._final_disk_read_bytes is not None
        ):
            disk_read_bytes = max(
                0,
                self._final_disk_read_bytes - self._initial_disk_read_bytes,
            )

        if (
            self._initial_disk_write_bytes is not None
            and self._final_disk_write_bytes is not None
        ):
            disk_write_bytes = max(
                0,
                self._final_disk_write_bytes - self._initial_disk_write_bytes,
            )

        return {
            "resource_metrics_available": self._process_handle is not None,
            "observed_max_rss_bytes": self._observed_max_rss_bytes,
            "disk_read_bytes": disk_read_bytes,
            "disk_write_bytes": disk_write_bytes,
        }

    def build_summary(self) -> Dict[str, Any]:
        """
        Build a JSON-serializable telemetry summary suitable for a manifest.
        """
        self.finalize()

        with self._state_lock:
            return {
                "telemetry_format_version": TELEMETRY_FORMAT_VERSION,
                "stages": {
                    stage_name: stage.build_summary()
                    for stage_name, stage in sorted(self._stages.items())
                },
                "local_phase_seconds": dict(
                    sorted(self._local_phase_seconds.items())
                ),
                "process": self._build_process_summary(),
            }


class NCBIRequestStartRateLimiter:
    """
    Enforce a maximum number of request *start* instants per second.

    NCBI's published allowance is expressed as requests per second, measured at
    the moment a request is issued rather than when it completes. A limiter that
    reserves start slots under one shared lock therefore bounds the aggregate
    rate correctly even when several worker threads are in flight, which a
    per-worker ``time.sleep`` cannot do.
    """

    def __init__(
        self,
        *,
        max_request_starts_per_second: float,
        sleep_function=time.sleep,
        monotonic_function=time.monotonic,
    ) -> None:
        if max_request_starts_per_second <= 0:
            raise ValueError(
                "max_request_starts_per_second must be a positive number."
            )

        self._state_lock = threading.Lock()
        self._sleep_function = sleep_function
        self._monotonic_function = monotonic_function
        self._minimum_start_interval_seconds = 1.0 / max_request_starts_per_second
        self._next_allowed_start_at_seconds = self._monotonic_function()
        self._max_request_starts_per_second = max_request_starts_per_second
        self.total_wait_seconds = 0.0
        self.acquired_start_count = 0

    @property
    def max_request_starts_per_second(self) -> float:
        with self._state_lock:
            return self._max_request_starts_per_second

    def set_max_request_starts_per_second(
        self,
        *,
        max_request_starts_per_second: float,
    ) -> None:
        if max_request_starts_per_second <= 0:
            raise ValueError(
                "max_request_starts_per_second must be a positive number."
            )

        with self._state_lock:
            self._max_request_starts_per_second = max_request_starts_per_second
            self._minimum_start_interval_seconds = (
                1.0 / max_request_starts_per_second
            )

    def acquire_request_start_slot(self) -> float:
        """
        Reserve the next permitted start instant and wait for it.

        Returns:
            float: seconds spent waiting for the reserved start slot.
        """
        with self._state_lock:
            current_time_seconds = self._monotonic_function()
            scheduled_start_at_seconds = max(
                current_time_seconds,
                self._next_allowed_start_at_seconds,
            )
            self._next_allowed_start_at_seconds = (
                scheduled_start_at_seconds + self._minimum_start_interval_seconds
            )
            self.acquired_start_count += 1

        wait_seconds = scheduled_start_at_seconds - self._monotonic_function()
        if wait_seconds > 0:
            self._sleep_function(wait_seconds)
            with self._state_lock:
                self.total_wait_seconds += wait_seconds
            return wait_seconds

        return 0.0


class NCBIAdaptiveConcurrencyGovernor:
    """
    Bound in-flight NCBI requests and back off when NCBI signals overload.

    The governor owns two coupled controls: how many requests may be in flight
    at once, and how fast new requests may start. Both are reduced together on a
    rate-limit, server-error, timeout, or circuit-breaker signal, and are
    recovered one step at a time only after a run of clean responses.
    """

    def __init__(
        self,
        *,
        max_concurrent_requests: int,
        rate_limiter: Optional[NCBIRequestStartRateLimiter] = None,
        minimum_concurrent_requests: int = 1,
        successes_before_recovery: int = 20,
        minimum_request_starts_per_second: float = 1.0,
    ) -> None:
        if max_concurrent_requests <= 0:
            raise ValueError("max_concurrent_requests must be a positive integer.")
        if minimum_concurrent_requests <= 0:
            raise ValueError(
                "minimum_concurrent_requests must be a positive integer."
            )
        if successes_before_recovery <= 0:
            raise ValueError(
                "successes_before_recovery must be a positive integer."
            )

        self._state_condition = threading.Condition(threading.Lock())
        self._configured_concurrent_requests = max_concurrent_requests
        self._minimum_concurrent_requests = min(
            minimum_concurrent_requests,
            max_concurrent_requests,
        )
        self._permitted_concurrent_requests = max_concurrent_requests
        self._in_flight_request_count = 0
        self._successes_before_recovery = successes_before_recovery
        self._consecutive_success_count = 0
        self._rate_limiter = rate_limiter
        self._minimum_request_starts_per_second = minimum_request_starts_per_second
        self._configured_request_starts_per_second = (
            rate_limiter.max_request_starts_per_second
            if rate_limiter is not None
            else None
        )
        self.reduction_event_count = 0
        self.recovery_event_count = 0

    @property
    def permitted_concurrent_requests(self) -> int:
        with self._state_condition:
            return self._permitted_concurrent_requests

    @contextmanager
    def request_slot(self):
        """
        Hold one in-flight request slot for the duration of a request attempt.
        """
        with self._state_condition:
            while self._in_flight_request_count >= self._permitted_concurrent_requests:
                self._state_condition.wait()
            self._in_flight_request_count += 1

        try:
            yield
        finally:
            with self._state_condition:
                self._in_flight_request_count -= 1
                self._state_condition.notify_all()

    def record_overload_signal(self) -> None:
        """
        Halve the permitted request rate and shrink concurrency by one step.
        """
        with self._state_condition:
            self._consecutive_success_count = 0
            previous_permitted_requests = self._permitted_concurrent_requests
            self._permitted_concurrent_requests = max(
                self._minimum_concurrent_requests,
                self._permitted_concurrent_requests - 1,
            )
            reduced = (
                self._permitted_concurrent_requests != previous_permitted_requests
            )
            self._state_condition.notify_all()

        if self._rate_limiter is not None:
            reduced_request_starts_per_second = max(
                self._minimum_request_starts_per_second,
                self._rate_limiter.max_request_starts_per_second / 2.0,
            )
            if (
                reduced_request_starts_per_second
                != self._rate_limiter.max_request_starts_per_second
            ):
                reduced = True
            self._rate_limiter.set_max_request_starts_per_second(
                max_request_starts_per_second=reduced_request_starts_per_second,
            )

        if reduced:
            with self._state_condition:
                self.reduction_event_count += 1

    def record_successful_request(self) -> None:
        """
        Count one clean response and restore one step once the run is long enough.
        """
        with self._state_condition:
            self._consecutive_success_count += 1
            if self._consecutive_success_count < self._successes_before_recovery:
                return

            self._consecutive_success_count = 0
            recovered = False
            if (
                self._permitted_concurrent_requests
                < self._configured_concurrent_requests
            ):
                self._permitted_concurrent_requests += 1
                self._state_condition.notify_all()
                recovered = True

        if (
            self._rate_limiter is not None
            and self._configured_request_starts_per_second is not None
        ):
            current_request_starts_per_second = (
                self._rate_limiter.max_request_starts_per_second
            )
            if (
                current_request_starts_per_second
                < self._configured_request_starts_per_second
            ):
                self._rate_limiter.set_max_request_starts_per_second(
                    max_request_starts_per_second=min(
                        self._configured_request_starts_per_second,
                        current_request_starts_per_second * 2.0,
                    ),
                )
                recovered = True

        if recovered:
            with self._state_condition:
                self.recovery_event_count += 1
