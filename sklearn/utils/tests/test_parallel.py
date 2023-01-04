import threading

import pytest
from joblib import Parallel

import sklearn
from sklearn.utils.fixes import delayed


def get_working_memory():
    return sklearn.get_config()["working_memory"]


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("backend", ["loky", "threading", "multiprocessing"])
def test_configuration_passes_through_to_joblib(n_jobs, backend):
    # Tests that the global global configuration is passed to joblib jobs

    n_iter, results = 10, []

    def parallel_inspect_config():
        with sklearn.config_context(working_memory=123):
            config = sklearn.get_config()
            results.extend(
                Parallel(n_jobs=n_jobs, pre_dispatch=2 * n_jobs, backend=backend)(
                    delayed(get_working_memory, config=config)() for _ in range(n_iter)
                )
            )

    other_thread = threading.Thread(target=parallel_inspect_config)
    other_thread.start()
    other_thread.join()

    assert results == [123] * n_iter
