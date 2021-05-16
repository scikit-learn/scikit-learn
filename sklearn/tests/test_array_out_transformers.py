import numpy as np
import pytest
from numpy.testing import assert_array_equal

from sklearn.utils._testing import ignore_warnings
from sklearn.base import clone
from sklearn.utils import all_estimators
from sklearn.utils._testing import set_random_state
from sklearn.utils._testing import SkipTest
from sklearn.utils._testing import raises
from sklearn.utils.estimator_checks import _enforce_estimator_tags_x
from sklearn.utils.estimator_checks import _enforce_estimator_tags_y
from sklearn.utils.estimator_checks import _pairwise_estimator_convert_X
from sklearn.utils.estimator_checks import _set_checking_parameters
from sklearn.utils.estimator_checks import _construct_instance
from sklearn.utils.estimator_checks import _get_check_estimator_ids
from sklearn.utils._tags import _safe_tags


def check_array_out_pandas(estimator_orig):
    """Check that array_out controls the output of a transformer."""
    try:
        import pandas as pd
    except ImportError:
        raise SkipTest("pandas is not installed: not testing for "
                       "array_out with pandas")

    tags = _safe_tags(estimator_orig)

    supported_tags = {"2darray", "categorical"}
    X_types_set = set(tags["X_types"])
    if not (supported_tags & X_types_set) or tags["no_validation"]:
        return

    estimator = clone(estimator_orig)
    set_random_state(estimator)
    _set_checking_parameters(estimator)

    if "warm_start" in estimator.get_params():
        estimator.set_params(warm_start=False)

    n_samples, n_features = 150, 8

    rng = np.random.RandomState(0)
    X = rng.normal(size=(n_samples, n_features))
    X = _enforce_estimator_tags_x(estimator, X)
    X = _pairwise_estimator_convert_X(X, estimator)

    if "categorical" in X_types_set:
        # simple transformation for categorical data
        X = (((X - X.min()) * 3) // 3).astype(np.int32)

    y = None
    if "requires_y" in tags:
        y = rng.randint(low=0, high=2, size=n_samples)
        y = _enforce_estimator_tags_y(estimator, y)

    feature_names = [f"feature_{i}" for i in range(X.shape[1])]
    X_train_df = pd.DataFrame(X, columns=feature_names)

    # Call `fit` on dataframe
    estimator.fit(X_train_df, y=y)
    assert all(estimator.feature_names_in_ == feature_names)

    has_transform = hasattr(estimator, "transform")
    has_fit_transform = hasattr(estimator, "fit_transform")

    test_index = [2*i for i in range(8)]
    X_test_df = pd.DataFrame(X[rng.permutation(8), :],
                             columns=feature_names, index=test_index)

    if has_transform:
        X_trans_df = estimator.transform(X_test_df, array_out='pandas')
        assert isinstance(X_trans_df, pd.DataFrame)
        assert_array_equal(X_trans_df.index, test_index)

    if has_fit_transform:
        X_trans_df = estimator.fit_transform(X_train_df, y=y,
                                             array_out='pandas')
        assert isinstance(X_trans_df, pd.DataFrame)
        assert_array_equal(X_trans_df.index, X_train_df.index)

    if has_transform:
        # Check that `transform` fails with dataframe with different names.
        X_test_invalid_df = pd.DataFrame(X[rng.permutation(8), :],
                                         columns=feature_names[::-1])
        match = ("The input's feature names does not match the "
                 "feature_names_in_ attribute.")
        with raises(ValueError, match=match):
            estimator.transform(X_test_invalid_df)

    # Fit on ndarray (without feature names)
    estimator.fit(X)

    if has_transform:
        X_trans_df = estimator.transform(X_test_df, array_out='pandas')
        assert isinstance(X_trans_df, pd.DataFrame)
        assert_array_equal(X_trans_df.index, test_index)

    if has_fit_transform:
        # Check that `fit_transform` also works
        X_trans_df = estimator.fit_transform(X, y=y, array_out='pandas')
        assert isinstance(X_trans_df, pd.DataFrame)
        assert_array_equal(X_trans_df.index, range(150))


ARRAY_OUT_TO_IGNORE = {
    'cluster',
    'compose',
    'cross_decomposition',
    'decomposition',
    'discriminant_analysis',
    'ensemble',
    'feature_extraction',
    'feature_selection',
    'impute',
    'isotonic',
    'kernel_approximation',
    'manifold',
    'neighbors',
    'neural_network',
    'pipeline',
    'random_projection',
}


def all_transformers_2d():
    with ignore_warnings(category=UserWarning):
        estimators = all_estimators(type_filter="transformer")

    for name, Estimator in estimators:
        module = Estimator.__module__.split(".")[1]
        if module in ARRAY_OUT_TO_IGNORE:
            continue

        try:
            estimator = _construct_instance(Estimator)
        except SkipTest:
            continue
        yield name, estimator


@ignore_warnings(category=UserWarning)
@pytest.mark.parametrize("name, estimator", all_transformers_2d(),
                         ids=_get_check_estimator_ids)
def test_array_out_pandas(name, estimator):
    check_array_out_pandas(estimator)
