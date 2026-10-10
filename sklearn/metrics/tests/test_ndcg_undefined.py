"""Tests for ndcg_score with replaced_undefined_by parameter (issue #29521)."""

import numpy as np
import pytest

from sklearn.exceptions import UndefinedMetricWarning
from sklearn.metrics import ndcg_score


def test_ndcg_score_replaced_undefined_by():
    """Check behavior of ndcg_score with replaced_undefined_by parameter."""
    y_true_undef = np.array([[0, 0, 0]])
    y_score = np.array([[0.1, 0.2, 0.3]])

    # Default: replaced_undefined_by=np.nan raises UndefinedMetricWarning
    with pytest.warns(UndefinedMetricWarning, match="NDCG is not defined"):
        result = ndcg_score(y_true_undef, y_score)
    assert np.isnan(result)

    # replaced_undefined_by=0.0: no warning, returns 0.0
    result = ndcg_score(y_true_undef, y_score, replaced_undefined_by=0.0)
    assert result == pytest.approx(0.0)

    # replaced_undefined_by=1.0: no warning, returns 1.0
    result = ndcg_score(y_true_undef, y_score, replaced_undefined_by=1.0)
    assert result == pytest.approx(1.0)

    # Mixed: one undefined sample, one defined
    y_true_mixed = np.array([[0, 0, 0], [1, 0, 1]])
    y_score_mixed = np.array([[0.1, 0.2, 0.3], [0.5, 0.1, 0.4]])

    # Default: warning raised, undefined sample excluded from average
    with pytest.warns(UndefinedMetricWarning, match="NDCG is not defined"):
        result_mixed = ndcg_score(y_true_mixed, y_score_mixed)
    result_defined_only = ndcg_score(y_true_mixed[1:], y_score_mixed[1:])
    assert result_mixed == pytest.approx(result_defined_only)

    # All samples undefined: returns nan
    y_true_all_zero = np.array([[0, 0, 0], [0, 0, 0]])
    with pytest.warns(UndefinedMetricWarning, match="NDCG is not defined"):
        result_all_undef = ndcg_score(y_true_all_zero, y_score_mixed)
    assert np.isnan(result_all_undef)

    # No undefined samples: no warning, finite result
    y_true_defined = np.array([[1, 0, 1], [0, 1, 0]])
    result_no_undef = ndcg_score(y_true_defined, y_score_mixed)
    assert np.isfinite(result_no_undef)


def test_ndcg_score_replaced_undefined_by_with_sample_weight():
    """Check replaced_undefined_by respects sample_weight."""
    # One undefined, one defined sample with explicit weights
    y_true = np.array([[0, 0, 0], [1, 0, 1]])
    y_score = np.array([[0.1, 0.2, 0.3], [0.5, 0.1, 0.4]])

    with pytest.warns(UndefinedMetricWarning):
        result_weighted = ndcg_score(y_true, y_score, sample_weight=[1.0, 2.0])
    result_unweighted = ndcg_score(y_true[1:], y_score[1:])
    # Both should give the NDCG of the one defined sample
    assert result_weighted == pytest.approx(result_unweighted)


def test_ndcg_score_replaced_undefined_by_validation():
    """Check replaced_undefined_by validates its input."""
    y_true = np.array([[1, 0, 1]])
    y_score = np.array([[0.5, 0.1, 0.4]])

    # Valid values: 0.0, 1.0, np.nan
    ndcg_score(y_true, y_score, replaced_undefined_by=0.0)
    ndcg_score(y_true, y_score, replaced_undefined_by=1.0)
    ndcg_score(y_true, y_score, replaced_undefined_by=np.nan)

    # Invalid values should raise
    with pytest.raises(ValueError):
        ndcg_score(y_true, y_score, replaced_undefined_by=1.5)
    with pytest.raises(ValueError):
        ndcg_score(y_true, y_score, replaced_undefined_by=-0.1)
