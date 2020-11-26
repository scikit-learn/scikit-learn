import pytest

from sklearn.base import BaseEstimator
from sklearn.utils._tags import (
    _DEFAULT_TAGS,
    _safe_tags,
)


class EstimatorNoTags:
    pass


class EstimatorOwnTags:
    def _get_tags(self):
        return {}


@pytest.mark.parametrize(
    "estimator, err_msg",
    [
        (BaseEstimator(), "The key xxx is neither defined"),
        (EstimatorNoTags(), "The key xxx is not a default tags defined"),
    ],
)
def test_safe_tags_error(estimator, err_msg):
    # Check that safe_tags raises error in ambiguous case.
    with pytest.raises(ValueError, match=err_msg):
        _safe_tags(estimator, key="xxx")


@pytest.mark.parametrize(
    "estimator, key, default, expected_results",
    [
        (EstimatorNoTags(), None, None, _DEFAULT_TAGS),
        (EstimatorNoTags(), "allow_nan", None, _DEFAULT_TAGS["allow_nan"]),
        (
            EstimatorNoTags(),
            "allow_nan",
            not _DEFAULT_TAGS["allow_nan"],
            not _DEFAULT_TAGS["allow_nan"],
        ),
        (EstimatorNoTags(), "xxx", True, True),
        (BaseEstimator(), None, None, _DEFAULT_TAGS),
        (BaseEstimator(), "allow_nan", None, _DEFAULT_TAGS["allow_nan"]),
        (
            BaseEstimator(),
            "allow_nan",
            not _DEFAULT_TAGS["allow_nan"],
            _DEFAULT_TAGS["allow_nan"],
        ),
        (BaseEstimator(), "xxx", True, True),
        (EstimatorOwnTags(), None, None, _DEFAULT_TAGS),
        (EstimatorOwnTags(), "allow_nan", None, _DEFAULT_TAGS["allow_nan"]),
        (
            EstimatorOwnTags(),
            "allow_nan",
            not _DEFAULT_TAGS["allow_nan"],
            not _DEFAULT_TAGS["allow_nan"],
        ),
        (EstimatorOwnTags(), "xxx", True, True),
    ],
)
def test_safe_tags_no_get_tags(estimator, key, default, expected_results):
    # check the behaviour of _safe_tags when an estimator does not implement
    # _get_tags
    assert _safe_tags(estimator, key=key, default=default) == expected_results
