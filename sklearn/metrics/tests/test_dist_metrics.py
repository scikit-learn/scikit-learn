import itertools
import pickle
import copy

import numpy as np
import pytest

import scipy.sparse as sp
from scipy.spatial.distance import cdist
from sklearn.metrics import DistanceMetric

from sklearn.metrics._dist_metrics import (
    BOOL_METRICS,
    # Unexposed private DistanceMetric for 32 bit
    DistanceMetric32,
)

from sklearn.utils import check_random_state
from sklearn.utils._testing import assert_allclose, create_memmap_backed_data
from sklearn.utils.fixes import sp_version, parse_version


def dist_func(x1, x2, p):
    return np.sum((x1 - x2) ** p) ** (1.0 / p)


rng = check_random_state(0)
d = 4
n1 = 20
n2 = 25
X64 = rng.random_sample((n1, d))
Y64 = rng.random_sample((n2, d))
X32 = X64.astype("float32")
Y32 = Y64.astype("float32")

[X_mmap, Y_mmap] = create_memmap_backed_data([X64, Y64])

# make boolean arrays: ones and zeros
X_bool = (X64 < 0.3).astype(np.float64)  # quite sparse
Y_bool = (Y64 < 0.7).astype(np.float64)  # not too sparse

[X_bool_mmap, Y_bool_mmap] = create_memmap_backed_data([X_bool, Y_bool])


V = rng.random_sample((d, d))
VI = np.dot(V, V.T)


METRICS_DEFAULT_PARAMS = [
    ("euclidean", {}),
    ("cityblock", {}),
    ("minkowski", dict(p=(1, 1.5, 2, 3))),
    ("chebyshev", {}),
    ("seuclidean", dict(V=(rng.random_sample(d),))),
    ("mahalanobis", dict(VI=(VI,))),
    ("hamming", {}),
    ("canberra", {}),
    ("braycurtis", {}),
]
if sp_version >= parse_version("1.8.0.dev0"):
    # Starting from scipy 1.8.0.dev0, minkowski now accepts w, the weighting
    # parameter directly and using it is preferred over using wminkowski.
    METRICS_DEFAULT_PARAMS.append(
        ("minkowski", dict(p=(1, 1.5, 3), w=(rng.random_sample(d),))),
    )
else:
    # For previous versions of scipy, this was possible through a dedicated
    # metric (deprecated in 1.6 and removed in 1.8).
    METRICS_DEFAULT_PARAMS.append(
        ("wminkowski", dict(p=(1, 1.5, 3), w=(rng.random_sample(d),))),
    )


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize(
    "metric_param_grid", METRICS_DEFAULT_PARAMS, ids=lambda params: params[0]
)
@pytest.mark.parametrize("X, Y", [(X64, Y64), (X32, Y32), (X_mmap, Y_mmap)])
def test_cdist(metric_param_grid, X, Y):
    DistanceMetricInterface = (
        DistanceMetric if X.dtype == Y.dtype == np.float64 else DistanceMetric32
    )
    metric, param_grid = metric_param_grid
    keys = param_grid.keys()
    X_csr, Y_csr = sp.csr_matrix(X), sp.csr_matrix(Y)
    for vals in itertools.product(*param_grid.values()):
        kwargs = dict(zip(keys, vals))
        rtol_dict = {}
        if metric == "mahalanobis" and X.dtype == np.float32:
            # Computation of mahalanobis differs between
            # the scipy and scikit-learn implementation.
            # Hence, we increase the relative tolerance.
            # TODO: Inspect slight numerical discrepancy
            # with scipy
            rtol_dict = {"rtol": 1e-6}

        if metric == "wminkowski":
            # wminkoski is deprecated in SciPy 1.6.0 and removed in 1.8.0
            WarningToExpect = None
            if sp_version >= parse_version("1.6.0"):
                WarningToExpect = DeprecationWarning
            with pytest.warns(WarningToExpect):
                D_scipy_cdist = cdist(X, Y, metric, **kwargs)
        else:
            D_scipy_cdist = cdist(X, Y, metric, **kwargs)

        dm = DistanceMetricInterface.get_metric(metric, **kwargs)

        # DistanceMetric.pairwise must be consistent for all
        # combinations of formats in {sparse, dense}.
        D_sklearn = dm.pairwise(X, Y)
        assert D_sklearn.flags.c_contiguous
        assert_allclose(D_sklearn, D_scipy_cdist, **rtol_dict)

        D_sklearn = dm.pairwise(X_csr, Y_csr)
        assert D_sklearn.flags.c_contiguous
        assert_allclose(D_sklearn, D_scipy_cdist, **rtol_dict)

        D_sklearn = dm.pairwise(X_csr, Y)
        assert D_sklearn.flags.c_contiguous
        assert_allclose(D_sklearn, D_scipy_cdist, **rtol_dict)

        D_sklearn = dm.pairwise(X, Y_csr)
        assert D_sklearn.flags.c_contiguous
        assert_allclose(D_sklearn, D_scipy_cdist, **rtol_dict)


@pytest.mark.parametrize("metric", BOOL_METRICS)
@pytest.mark.parametrize(
    "X_bool, Y_bool", [(X_bool, Y_bool), (X_bool_mmap, Y_bool_mmap)]
)
def test_cdist_bool_metric(metric, X_bool, Y_bool):
    D_scipy_cdist = cdist(X_bool, Y_bool, metric)

    dm = DistanceMetric.get_metric(metric)
    D_sklearn = dm.pairwise(X_bool, Y_bool)
    assert_allclose(D_sklearn, D_scipy_cdist)

    # DistanceMetric.pairwise must be consistent
    # on all combinations of format in {sparse, dense}².
    X_bool_csr, Y_bool_csr = sp.csr_matrix(X_bool), sp.csr_matrix(Y_bool)

    D_sklearn = dm.pairwise(X_bool, Y_bool)
    assert D_sklearn.flags.c_contiguous
    assert_allclose(D_sklearn, D_scipy_cdist)

    D_sklearn = dm.pairwise(X_bool_csr, Y_bool_csr)
    assert D_sklearn.flags.c_contiguous
    assert_allclose(D_sklearn, D_scipy_cdist)

    D_sklearn = dm.pairwise(X_bool, Y_bool_csr)
    assert D_sklearn.flags.c_contiguous
    assert_allclose(D_sklearn, D_scipy_cdist)

    D_sklearn = dm.pairwise(X_bool_csr, Y_bool)
    assert D_sklearn.flags.c_contiguous
    assert_allclose(D_sklearn, D_scipy_cdist)


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize(
    "metric_param_grid", METRICS_DEFAULT_PARAMS, ids=lambda params: params[0]
)
@pytest.mark.parametrize("X", [X64, X32, X_mmap])
def test_pdist(metric_param_grid, X):
    DistanceMetricInterface = (
        DistanceMetric if X.dtype == np.float64 else DistanceMetric32
    )
    metric, param_grid = metric_param_grid
    keys = param_grid.keys()
    X_csr = sp.csr_matrix(X)
    for vals in itertools.product(*param_grid.values()):
        kwargs = dict(zip(keys, vals))
        rtol_dict = {}
        if metric == "mahalanobis" and X.dtype == np.float32:
            # Computation of mahalanobis differs between
            # the scipy and scikit-learn implementation.
            # Hence, we increase the relative tolerance.
            # TODO: Inspect slight numerical discrepancy
            # with scipy
            rtol_dict = {"rtol": 1e-6}

        if metric == "wminkowski":
            if sp_version >= parse_version("1.8.0"):
                pytest.skip("wminkowski will be removed in SciPy 1.8.0")

            # wminkoski is deprecated in SciPy 1.6.0 and removed in 1.8.0
            ExceptionToAssert = None
            if sp_version >= parse_version("1.6.0"):
                ExceptionToAssert = DeprecationWarning
            with pytest.warns(ExceptionToAssert):
                D_scipy_pdist = cdist(X, X, metric, **kwargs)
        else:
            D_scipy_pdist = cdist(X, X, metric, **kwargs)

        dm = DistanceMetricInterface.get_metric(metric, **kwargs)
        D_sklearn = dm.pairwise(X)
        assert D_sklearn.flags.c_contiguous
        assert_allclose(D_sklearn, D_scipy_pdist, **rtol_dict)

        D_sklearn_csr = dm.pairwise(X_csr)
        assert D_sklearn.flags.c_contiguous
        assert_allclose(D_sklearn_csr, D_scipy_pdist, **rtol_dict)

        D_sklearn_csr = dm.pairwise(X_csr, X_csr)
        assert D_sklearn.flags.c_contiguous
        assert_allclose(D_sklearn_csr, D_scipy_pdist, **rtol_dict)


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize(
    "metric_param_grid", METRICS_DEFAULT_PARAMS, ids=lambda params: params[0]
)
def test_distance_metrics_dtype_consistency(metric_param_grid):
    # DistanceMetric must return similar distances for both float32 and float64
    # input data.
    metric, param_grid = metric_param_grid
    keys = param_grid.keys()

    # Choose rtol to make sure that this test is robust to changes in the random
    # seed in the module-level test data generation code.
    rtol = 1e-5

    for vals in itertools.product(*param_grid.values()):
        kwargs = dict(zip(keys, vals))
        dm64 = DistanceMetric.get_metric(metric, **kwargs)
        dm32 = DistanceMetric32.get_metric(metric, **kwargs)

        D64 = dm64.pairwise(X64)
        D32 = dm32.pairwise(X32)

        # Both results are np.float64 dtype because the accumulation accross
        # features is done in float64. However the input data and the element
        # wise arithmetic operations are done in float32 so we can expect a
        # small discrepancy.
        assert D64.dtype == D32.dtype == np.float64

        # assert_allclose introspects the dtype of the input arrays to decide
        # which rtol value to use by default but in this case we know that D32
        # is not computed with the same precision so we set rtol manually.
        assert_allclose(D64, D32, rtol=rtol)

        D64 = dm64.pairwise(X64, Y64)
        D32 = dm32.pairwise(X32, Y32)
        assert_allclose(D64, D32, rtol=rtol)


@pytest.mark.parametrize("metric", BOOL_METRICS)
@pytest.mark.parametrize("X_bool", [X_bool, X_bool_mmap])
def test_pdist_bool_metrics(metric, X_bool):
    D_scipy_pdist = cdist(X_bool, X_bool, metric)
    dm = DistanceMetric.get_metric(metric)
    D_sklearn = dm.pairwise(X_bool)
    assert_allclose(D_sklearn, D_scipy_pdist)

    X_bool_csr = sp.csr_matrix(X_bool)
    D_sklearn = dm.pairwise(X_bool_csr)
    assert_allclose(D_sklearn, D_scipy_pdist)


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize("writable_kwargs", [True, False])
@pytest.mark.parametrize(
    "metric_param_grid", METRICS_DEFAULT_PARAMS, ids=lambda params: params[0]
)
@pytest.mark.parametrize("X", [X64, X32])
def test_pickle(writable_kwargs, metric_param_grid, X):
    DistanceMetricInterface = (
        DistanceMetric if X.dtype == np.float64 else DistanceMetric32
    )
    metric, param_grid = metric_param_grid
    keys = param_grid.keys()
    for vals in itertools.product(*param_grid.values()):
        if any(isinstance(val, np.ndarray) for val in vals):
            vals = copy.deepcopy(vals)
            for val in vals:
                if isinstance(val, np.ndarray):
                    val.setflags(write=writable_kwargs)
        kwargs = dict(zip(keys, vals))
        dm = DistanceMetricInterface.get_metric(metric, **kwargs)
        D1 = dm.pairwise(X)
        dm2 = pickle.loads(pickle.dumps(dm))
        D2 = dm2.pairwise(X)
        assert_allclose(D1, D2)


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize("metric", BOOL_METRICS)
@pytest.mark.parametrize("X_bool", [X_bool, X_bool_mmap])
def test_pickle_bool_metrics(metric, X_bool):
    dm = DistanceMetric.get_metric(metric)
    D1 = dm.pairwise(X_bool)
    dm2 = pickle.loads(pickle.dumps(dm))
    D2 = dm2.pairwise(X_bool)
    assert_allclose(D1, D2)


@pytest.mark.parametrize("X, Y", [(X64, Y64), (X32, Y32), (X_mmap, Y_mmap)])
def test_haversine_metric(X, Y):
    DistanceMetricInterface = (
        DistanceMetric if X.dtype == np.float64 else DistanceMetric32
    )

    # The Haversine DistanceMetric only works on 2 features.
    X = np.asarray(X[:, :2])
    Y = np.asarray(Y[:, :2])

    X_csr, Y_csr = sp.csr_matrix(X), sp.csr_matrix(Y)

    # Haversine is not supported by scipy.special.distance.{cdist,pdist}
    # So we reimplement it to have a reference.
    def haversine_slow(x1, x2):
        return 2 * np.arcsin(
            np.sqrt(
                np.sin(0.5 * (x1[0] - x2[0])) ** 2
                + np.cos(x1[0]) * np.cos(x2[0]) * np.sin(0.5 * (x1[1] - x2[1])) ** 2
            )
        )

    D_reference = np.zeros((X_csr.shape[0], Y_csr.shape[0]))
    for i, xi in enumerate(X):
        for j, yj in enumerate(Y):
            D_reference[i, j] = haversine_slow(xi, yj)

    haversine = DistanceMetricInterface.get_metric("haversine")

    D_sklearn = haversine.pairwise(X, Y)
    assert_allclose(
        haversine.dist_to_rdist(D_sklearn), np.sin(0.5 * D_reference) ** 2, rtol=1e-6
    )

    assert_allclose(D_sklearn, D_reference)

    D_sklearn = haversine.pairwise(X_csr, Y_csr)
    assert D_sklearn.flags.c_contiguous
    assert_allclose(D_sklearn, D_reference)

    D_sklearn = haversine.pairwise(X_csr, Y)
    assert D_sklearn.flags.c_contiguous
    assert_allclose(D_sklearn, D_reference)

    D_sklearn = haversine.pairwise(X, Y_csr)
    assert D_sklearn.flags.c_contiguous
    assert_allclose(D_sklearn, D_reference)


def test_pyfunc_metric():
    X = np.random.random((10, 3))

    euclidean = DistanceMetric.get_metric("euclidean")
    pyfunc = DistanceMetric.get_metric("pyfunc", func=dist_func, p=2)

    # Check if both callable metric and predefined metric initialized
    # DistanceMetric object is picklable
    euclidean_pkl = pickle.loads(pickle.dumps(euclidean))
    pyfunc_pkl = pickle.loads(pickle.dumps(pyfunc))

    D1 = euclidean.pairwise(X)
    D2 = pyfunc.pairwise(X)

    D1_pkl = euclidean_pkl.pairwise(X)
    D2_pkl = pyfunc_pkl.pairwise(X)

    assert_allclose(D1, D2)
    assert_allclose(D1_pkl, D2_pkl)


def test_input_data_size():
    # Regression test for #6288
    # Previously, a metric requiring a particular input dimension would fail
    def custom_metric(x, y):
        assert x.shape[0] == 3
        return np.sum((x - y) ** 2)

    rng = check_random_state(0)
    X = rng.rand(10, 3)

    pyfunc = DistanceMetric.get_metric("pyfunc", func=custom_metric)
    eucl = DistanceMetric.get_metric("euclidean")
    assert_allclose(pyfunc.pairwise(X), eucl.pairwise(X) ** 2)


# TODO: Remove filterwarnings in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
def test_readonly_kwargs():
    # Non-regression test for:
    # https://github.com/scikit-learn/scikit-learn/issues/21685

    rng = check_random_state(0)

    weights = rng.rand(100)
    VI = rng.rand(10, 10)
    weights.setflags(write=False)
    VI.setflags(write=False)

    # Those distances metrics have to support readonly buffers.
    DistanceMetric.get_metric("seuclidean", V=weights)
    DistanceMetric.get_metric("wminkowski", p=1, w=weights)
    DistanceMetric.get_metric("mahalanobis", VI=VI)


@pytest.mark.parametrize(
    "w, err_type, err_msg",
    [
        (np.array([1, 1.5, -13]), ValueError, "w cannot contain negative weights"),
        (np.array([1, 1.5, np.nan]), ValueError, "w contains NaN"),
        (
            sp.csr_matrix([1, 1.5, 1]),
            TypeError,
            "A sparse matrix was passed, but dense data is required",
        ),
        (np.array(["a", "b", "c"]), ValueError, "could not convert string to float"),
        (np.array([]), ValueError, "a minimum of 1 is required"),
    ],
)
def test_minkowski_metric_validate_weights_values(w, err_type, err_msg):
    with pytest.raises(err_type, match=err_msg):
        DistanceMetric.get_metric("minkowski", p=3, w=w)


def test_minkowski_metric_validate_weights_size():
    w2 = rng.random_sample(d + 1)
    dm = DistanceMetric.get_metric("minkowski", p=3, w=w2)
    msg = (
        "MinkowskiDistance: the size of w must match "
        f"the number of features \\({X64.shape[1]}\\). "
        f"Currently len\\(w\\)={w2.shape[0]}."
    )
    with pytest.raises(ValueError, match=msg):
        dm.pairwise(X64, Y64)


# TODO: Remove in 1.3 when wminkowski is removed
def test_wminkowski_deprecated():
    w = rng.random_sample(d)
    msg = "WMinkowskiDistance is deprecated in version 1.1"
    with pytest.warns(FutureWarning, match=msg):
        DistanceMetric.get_metric("wminkowski", p=3, w=w)


# TODO: Remove in 1.3 when wminkowski is removed
@pytest.mark.filterwarnings("ignore:WMinkowskiDistance:FutureWarning:sklearn")
@pytest.mark.parametrize("p", [1, 1.5, 3])
def test_wminkowski_minkowski_equivalence(p):
    w = rng.random_sample(d)
    # Weights are rescaled for consistency w.r.t scipy 1.8 refactoring of 'minkowski'
    dm_wmks = DistanceMetric.get_metric("wminkowski", p=p, w=(w) ** (1 / p))
    dm_mks = DistanceMetric.get_metric("minkowski", p=p, w=w)
    D_wmks = dm_wmks.pairwise(X64, Y64)
    D_mks = dm_mks.pairwise(X64, Y64)
    assert_allclose(D_wmks, D_mks)
