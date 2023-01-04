# Authors: Gael Varoquaux <gael.varoquaux@normalesup.org>
#          Justin Vincent
#          Lars Buitinck
# License: BSD 3 clause

import math
import time

import numpy as np
import pytest
import scipy.stats
from joblib import Parallel

from sklearn import config_context
from sklearn.compose import make_column_transformer
from sklearn.datasets import load_iris
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConfigPropagationWarning
from sklearn.model_selection import GridSearchCV
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.utils._testing import assert_array_equal

from sklearn.utils.fixes import _object_dtype_isnan, delayed, loguniform


@pytest.mark.parametrize("dtype, val", ([object, 1], [object, "a"], [float, 1]))
def test_object_dtype_isnan(dtype, val):
    X = np.array([[val, np.nan], [np.nan, val]], dtype=dtype)

    expected_mask = np.array([[False, True], [True, False]])

    mask = _object_dtype_isnan(X)

    assert_array_equal(mask, expected_mask)


@pytest.mark.parametrize("low,high,base", [(-1, 0, 10), (0, 2, np.exp(1)), (-1, 1, 2)])
def test_loguniform(low, high, base):
    rv = loguniform(base**low, base**high)
    assert isinstance(rv, scipy.stats._distn_infrastructure.rv_frozen)
    rvs = rv.rvs(size=2000, random_state=0)

    # Test the basics; right bounds, right size
    assert (base**low <= rvs).all() and (rvs <= base**high).all()
    assert len(rvs) == 2000

    # Test that it's actually (fairly) uniform
    log_rvs = np.array([math.log(x, base) for x in rvs])
    counts, _ = np.histogram(log_rvs)
    assert counts.mean() == 200
    assert np.abs(counts - counts.mean()).max() <= 40

    # Test that random_state works
    assert loguniform(base**low, base**high).rvs(random_state=0) == loguniform(
        base**low, base**high
    ).rvs(random_state=0)


def test_delayed_warning_config():
    """Check that we raise an informing warning when the user does not pass
    the configuration to be used within the threads executing the tasks.
    """

    n_tasks = 10
    warning_msg = "scikit-learn 1.3 and above will require a config parameter to be"
    with pytest.warns(ConfigPropagationWarning, match=warning_msg) as records:
        Parallel(n_jobs=2, pre_dispatch=1)(
            delayed(time.sleep)(0) for _ in range(n_tasks)
        )

    assert len(records) == n_tasks


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
