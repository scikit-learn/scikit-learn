# import statements
import numpy as np
import pytest

# sklearn imports
from sklearn.utils._multiindexset import (
    FullTensor,
    HyperbolicCross,
    MultiIndexSet,
    TotalDegree,
    ZarembaCross,
)


# check multiindex set sizes
@pytest.mark.parametrize(
    "multiindex_set, size",
    [
        (FullTensor, 512),
        (TotalDegree, 120),
        (HyperbolicCross, 38),
        (ZarembaCross, 98),
    ],
)
def test_indices(multiindex_set, size):
    # dimension 3, degree 7
    multiindex_set = multiindex_set(3, 7)

    # generate indices
    assert len(list(multiindex_set.indices())) == size


# check weighted multiindex set sizes
@pytest.mark.parametrize(
    "multiindex_set, size",
    [
        (FullTensor, 28),
        (TotalDegree, 17),
        (HyperbolicCross, 11),
        (ZarembaCross, 15),
    ],
)
def test_indices_weighted(multiindex_set, size):
    # dimension 2, degree 6, weights (1, 3/5)
    multiindex_set = multiindex_set(2, 6, weights=(1, 3 / 5))

    # generate indices
    assert len(list(multiindex_set.indices())) == size


# check dimension
def test_dimension():
    # non-int dimension throws error
    with pytest.raises(ValueError, match="wololo"):
        TotalDegree("wololo", 2)

    # dimension < 0 throws error
    with pytest.raises(ValueError, match="-1"):
        TotalDegree(-1, 2)

    # dimension = 0 throws error
    with pytest.raises(ValueError, match="0"):
        TotalDegree(0, 2)

    # passes
    TotalDegree(1, 2)


# check degree
def test_degree():
    # non-int degree throws error
    with pytest.raises(ValueError, match="wololo"):
        TotalDegree(3, "wololo")

    # degree < 0 throws error
    with pytest.raises(ValueError, match="-1"):
        TotalDegree(3, -1)

    # degree = 0 passes
    TotalDegree(3, 0)

    # passes
    TotalDegree(3, 1)


# checks on weights
def test_weights():
    # weights with different type throws error
    with pytest.raises(ValueError, match="weights"):
        TotalDegree(3, 2, weights=1)

    # weights with different len throws error
    with pytest.raises(ValueError, match="weights"):
        TotalDegree(3, 2, weights=(0.5, 1))

    # weights that do not supper len throw error
    with pytest.raises(ValueError, match="weights"):
        TotalDegree(3, 2, weights="abc")

    # non-numeric weights throw error
    with pytest.raises(ValueError, match="weights"):
        TotalDegree(3, 2, weights=(0.5, 1, "wololo"))

    # weights must be > 0
    with pytest.raises(ValueError, match="weights"):
        TotalDegree(3, 2, weights=(0.5, 1, 0))

    # 2d array throws error
    with pytest.raises(ValueError, match="weights"):
        TotalDegree(4, 2, weights=np.array([[1, 2], [3, 4]]))

    # weights as 1d array passes
    TotalDegree(3, 2, weights=np.array([1, 2, 3]))

    # weights as generator passes
    TotalDegree(3, 2, weights=(j for j in range(1, 4)))


# special checks for Zaremba cross
def test_Zaremba():
    # test degree 0
    len(list(ZarembaCross(3, 0).indices())) == 1


# test lookup from string
@pytest.mark.parametrize(
    "name, upper_case_name, multiindex_set",
    [
        ("full_tensor", "FullTensor", FullTensor),
        ("total_degree", "TotalDegree", TotalDegree),
        ("hyperbolic_cross", "HyperbolicCross", HyperbolicCross),
        ("Zaremba_cross", "ZarembaCross", ZarembaCross),
    ],
)
def test_from_string(name, upper_case_name, multiindex_set):
    # unknown multiindex set type throws error
    with pytest.raises(ValueError, match="type"):
        MultiIndexSet.from_string("wololo")

    # UpperCase multiindex set throws error
    with pytest.raises(ValueError, match="type"):
        MultiIndexSet.from_string(upper_case_name)

    # passes
    assert MultiIndexSet.from_string(name) == multiindex_set


# test print method
def test_print():
    # unweighted
    multiindex_set = TotalDegree(2, 3)
    assert "3" in str(multiindex_set)
    assert "2" in str(multiindex_set)

    # weighted
    multiindex_set = TotalDegree(2, 3, weights=(1, 3 / 5))
    assert "weights" in str(multiindex_set)
