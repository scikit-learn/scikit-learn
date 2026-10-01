import os

import pytest

from sklearn.utils._openmp_helpers import (
    _openmp_effective_n_threads,
    _openmp_parallelism_enabled,
)


@pytest.mark.thread_unsafe  # changes schedaffinity, which is process global
@pytest.mark.skipif(
    not hasattr(os, "sched_setaffinity") or not _openmp_parallelism_enabled(),
    reason="Missing OS function and/or OpenMP",
)
@pytest.mark.skipif(
    _openmp_effective_n_threads() == 1, reason="Need access to more than one core"
)
@pytest.mark.skipif(
    os.sched_getaffinity(0) == 1, reason="Need access to more than one core"
)
@pytest.mark.skipif(
    "OMP_NUM_THREADS" in os.environ, reason="Overriding env variable set"
)
def test_openmp_cpus_available_notices_changes(request):
    """When available cores change, ``_openmp_effective_n_threads()`` notices."""
    current_affinity = os.sched_getaffinity(0)
    request.addfinalizer(lambda: os.sched_setaffinity(0, current_affinity))

    os.sched_setaffinity(0, list(current_affinity)[:1])
    assert _openmp_effective_n_threads() == 1
