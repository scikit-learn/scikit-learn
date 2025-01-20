import time
import warnings

import joblib
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn import config_context, get_config
from sklearn.compose import make_column_transformer
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils.fixes import _IS_WASM
from sklearn.utils.parallel import Parallel, delayed


def get_working_memory():
    return get_config()["working_memory"]


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("backend", ["loky", "threading", "multiprocessing"])
def test_configuration_passes_through_to_joblib(n_jobs, backend):
    # Tests that the global global configuration is passed to joblib jobs

    with config_context(working_memory=123):
        results = Parallel(n_jobs=n_jobs, backend=backend)(
            delayed(get_working_memory)() for _ in range(2)
        )

    assert_array_equal(results, [123] * 2)


def test_parallel_delayed_warnings():
    """Informative warnings should be raised when mixing sklearn and joblib API"""
    # We should issue a warning when one wants to use sklearn.utils.fixes.Parallel
    # with joblib.delayed. The config will not be propagated to the workers.
    warn_msg = "`sklearn.utils.parallel.Parallel` needs to be used in conjunction"
    with pytest.warns(UserWarning, match=warn_msg) as records:
        Parallel()(joblib.delayed(time.sleep)(0) for _ in range(10))
    assert len(records) == 10

    # We should issue a warning if one wants to use sklearn.utils.fixes.delayed with
    # joblib.Parallel
    warn_msg = (
        "`sklearn.utils.parallel.delayed` should be used with "
        "`sklearn.utils.parallel.Parallel` to make it possible to propagate"
    )
    with pytest.warns(UserWarning, match=warn_msg) as records:
        joblib.Parallel()(delayed(time.sleep)(0) for _ in range(10))
    assert len(records) == 10


@pytest.mark.parametrize("n_jobs", [1, 2])
def test_dispatch_config_parallel(n_jobs):
    """Check that we properly dispatch the configuration in parallel processing.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/25239
    """
    pd = pytest.importorskip("pandas")
    iris = load_iris(as_frame=True)

    class TransformerRequiredDataFrame(StandardScaler):
        def fit(self, X, y=None):
            assert isinstance(X, pd.DataFrame), "X should be a DataFrame"
            return super().fit(X, y)

        def transform(self, X, y=None):
            assert isinstance(X, pd.DataFrame), "X should be a DataFrame"
            return super().transform(X, y)

    dropper = make_column_transformer(
        ("drop", [0]),
        remainder="passthrough",
        n_jobs=n_jobs,
    )
    param_grid = {"randomforestclassifier__max_depth": [1, 2, 3]}
    search_cv = GridSearchCV(
        make_pipeline(
            dropper,
            TransformerRequiredDataFrame(),
            RandomForestClassifier(n_estimators=5, n_jobs=n_jobs),
        ),
        param_grid,
        cv=5,
        n_jobs=n_jobs,
        error_score="raise",  # this search should not fail
    )

    # make sure that `fit` would fail in case we don't request dataframe
    with pytest.raises(AssertionError, match="X should be a DataFrame"):
        search_cv.fit(iris.data, iris.target)

    with config_context(transform_output="pandas"):
        # we expect each intermediate steps to output a DataFrame
        search_cv.fit(iris.data, iris.target)

    assert not np.isnan(search_cv.cv_results_["mean_test_score"]).any()


def raise_warning():
    warnings.warn("Convergence warning", ConvergenceWarning)


@pytest.mark.parametrize("n_jobs", [1, 2])
@pytest.mark.parametrize("backend", ["loky", "threading", "multiprocessing"])
def test_filter_warning_propagates(n_jobs, backend):
    """Check warning propagates to the job."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", category=ConvergenceWarning)

        with pytest.raises(ConvergenceWarning):
            Parallel(n_jobs=n_jobs, backend=backend)(
                delayed(raise_warning)() for _ in range(2)
            )


def get_warnings():
    return warnings.filters


def test_check_warnings_threading():
    """Check that warnings filters are set correctly in the threading backend."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", category=ConvergenceWarning)

        filters = warnings.filters
        assert ("error", None, ConvergenceWarning, None, 0) in filters

        all_warnings = Parallel(n_jobs=2, backend="threading")(
            delayed(get_warnings)() for _ in range(2)
        )

        assert all(w == filters for w in all_warnings)


@pytest.mark.xfail(_IS_WASM, reason="Pyodide always use the sequential backend")
def test_filter_warning_propagates_no_side_effect_with_loky_backend():
    with warnings.catch_warnings():
        warnings.simplefilter("error", category=ConvergenceWarning)

        Parallel(n_jobs=2, backend="loky")(delayed(time.sleep)(0) for _ in range(10))

        # Since loky workers are reused, make sure that inside the loky workers,
        # warnings filters have been reset to their original value. Using joblib
        # directly should not turn ConvergenceWarning into an error.
        joblib.Parallel(n_jobs=2, backend="loky")(
            joblib.delayed(warnings.warn)("Convergence warning", ConvergenceWarning)
            for _ in range(10)
        )
