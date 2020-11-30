import pytest

from sklearn.base import BaseEstimator
from sklearn.utils._tags import (
    _DEFAULT_TAGS,
    _safe_tags,
)


class NoTagsEstimator:
    pass


@pytest.mark.parametrize(
    "estimator, err_msg",
    [
        (BaseEstimator(), "The key xxx is neither defined"),
        (NoTagsEstimator(), "The key xxx is not a default tags defined"),
    ],
)
def test_safe_tags_error(estimator, err_msg):
    # Check that safe_tags raises error in ambiguous case.
    with pytest.raises(ValueError, match=err_msg):
        _safe_tags(estimator, key="xxx")


@pytest.mark.parametrize(
    "estimator, key, default, expected_results",
    [
        (NoTagsEstimator(), None, None, _DEFAULT_TAGS),
        (NoTagsEstimator(), "allow_nan", None, _DEFAULT_TAGS["allow_nan"]),
        (
            NoTagsEstimator(),
            "allow_nan",
            not _DEFAULT_TAGS["allow_nan"],
            not _DEFAULT_TAGS["allow_nan"],
        ),
        (NoTagsEstimator(), "xxx", True, True),
        (BaseEstimator(), None, None, _DEFAULT_TAGS),
        (BaseEstimator(), "allow_nan", None, _DEFAULT_TAGS["allow_nan"]),
        (
            BaseEstimator(),
            "allow_nan",
            not _DEFAULT_TAGS["allow_nan"],
            _DEFAULT_TAGS["allow_nan"],
        ),
    ],
)
def test_safe_tags_no_get_tags(estimator, key, default, expected_results):
    # check the behaviour of _safe_tags when an estimator does not implement
    # _get_tags
    assert _safe_tags(estimator, key=key, default=default) == expected_results
