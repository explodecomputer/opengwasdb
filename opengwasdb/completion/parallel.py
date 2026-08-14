"""Process-pool worker setup shared by dense and ragged block-parallel completion."""
from __future__ import annotations

from typing import Any

from threadpoolctl import threadpool_limits  # type: ignore[import-untyped]

# Kept alive for the worker's lifetime so the BLAS thread cap persists (a bare
# threadpool_limits() call would reset on garbage collection).
_worker_thread_limiter: Any = None


def init_block_worker() -> None:
    """Cap each pool worker's BLAS (OpenBLAS/MKL) to one thread. numpy's linear
    algebra and sklearn otherwise spawn one thread per core *inside every worker*,
    so n_workers processes each with ~n_core threads massively oversubscribe the
    CPU (~1000 threads on 256 cores). One BLAS thread per worker gives clean
    process-level parallelism — scale with n_workers, not threads."""
    global _worker_thread_limiter
    _worker_thread_limiter = threadpool_limits(limits=1)
