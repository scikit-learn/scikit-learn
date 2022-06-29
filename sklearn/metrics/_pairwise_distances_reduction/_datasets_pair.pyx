import numpy as np
cimport numpy as cnp

from cython cimport final
from scipy.sparse import issparse

from ...utils._typedefs cimport DTYPE_t, ITYPE_t
from ...metrics._dist_metrics cimport DistanceMetric

cnp.import_array()

cdef class DatasetsPair:
    """Abstract class which wraps a pair of datasets (X, Y).

    This class allows computing distances between a single pair of rows of
    of X and Y at a time given the pair of their indices (i, j). This class is
    specialized for each metric thanks to the :func:`get_for` factory classmethod.

    The handling of parallelization over chunks to compute the distances
    and aggregation for several rows at a time is done in dedicated
    subclasses of PairwiseDistancesReduction that in-turn rely on
    subclasses of DatasetsPair for each pair of rows in the data. The goal
    is to make it possible to decouple the generic parallelization and
    aggregation logic from metric-specific computation as much as
    possible.

    X and Y can be stored as C-contiguous np.ndarrays or CSR matrices
    in subclasses.

    This class avoids the overhead of dispatching distance computations
    to :class:`sklearn.metrics.DistanceMetric` based on the physical
    representation of the vectors (sparse vs. dense). It makes use of
    cython.final to remove the overhead of dispatching method calls.

    Parameters
    ----------
    distance_metric: DistanceMetric
        The distance metric responsible for computing distances
        between two vectors of (X, Y).
    """

    @classmethod
    def get_for(
        cls,
        X,
        Y,
        str metric="euclidean",
        dict metric_kwargs=None,
    ) -> DatasetsPair:
        """Return the DatasetsPair implementation for the given arguments.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples_X, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.
            If provided as a sparse matrix, it must be in CSR format.

        Y : {ndarray, sparse matrix} of shape (n_samples_Y, n_features)
            Input data.
            If provided as a ndarray, it must be C-contiguous.
            If provided as a sparse matrix, it must be in CSR format.

        metric : str, default='euclidean'
            The distance metric to compute between rows of X and Y.
            The default metric is a fast implementation of the Euclidean
            metric. For a list of available metrics, see the documentation
            of :class:`~sklearn.metrics.DistanceMetric`.

        metric_kwargs : dict, default=None
            Keyword arguments to pass to specified metric function.

        Returns
        -------
        datasets_pair: DatasetsPair
            The suited DatasetsPair implementation.
        """
        cdef:
            DistanceMetric distance_metric = DistanceMetric.get_metric(
                metric,
                **(metric_kwargs or {})
            )

        if not(X.dtype == Y.dtype == np.float64):
            raise ValueError(
                f"Only 64bit float datasets are supported at this time, "
                f"got: X.dtype={X.dtype} and Y.dtype={Y.dtype}."
            )

        # Metric-specific checks that do not replace nor duplicate `check_array`.
        distance_metric._validate_data(X)
        distance_metric._validate_data(Y)

        # TODO: dispatch to other dataset pairs for sparse support once available:
        if issparse(X) or issparse(Y):
            raise ValueError("Only dense datasets are supported for X and Y.")

        return DenseDenseDatasetsPair(X, Y, distance_metric)

    def __init__(self, DistanceMetric distance_metric):
        self.distance_metric = distance_metric

    cdef ITYPE_t n_samples_X(self) nogil:
        """Number of samples in X."""
        # This is a abstract method.
        # This _must_ always be overwritten in subclasses.
        # TODO: add "with gil: raise" here when supporting Cython 3.0
        return -999

    cdef ITYPE_t n_samples_Y(self) nogil:
        """Number of samples in Y."""
        # This is a abstract method.
        # This _must_ always be overwritten in subclasses.
        # TODO: add "with gil: raise" here when supporting Cython 3.0
        return -999

    cdef DTYPE_t surrogate_dist(self, ITYPE_t i, ITYPE_t j) nogil:
        return self.dist(i, j)

    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil:
        # This is a abstract method.
        # This _must_ always be overwritten in subclasses.
        # TODO: add "with gil: raise" here when supporting Cython 3.0
        return -1

@final
cdef class DenseDenseDatasetsPair(DatasetsPair):
    """Compute distances between row vectors of two arrays.

    Parameters
    ----------
    X: ndarray of shape (n_samples_X, n_features)
        Rows represent vectors. Must be C-contiguous.

    Y: ndarray of shape (n_samples_Y, n_features)
        Rows represent vectors. Must be C-contiguous.

    distance_metric: DistanceMetric
        The distance metric responsible for computing distances
        between two row vectors of (X, Y).
    """

    def __init__(self, X, Y, DistanceMetric distance_metric):
        super().__init__(distance_metric)
        # Arrays have already been checked
        self.X = X
        self.Y = Y
        self.d = X.shape[1]

    @final
    cdef ITYPE_t n_samples_X(self) nogil:
        return self.X.shape[0]

    @final
    cdef ITYPE_t n_samples_Y(self) nogil:
        return self.Y.shape[0]

    @final
    cdef DTYPE_t surrogate_dist(self, ITYPE_t i, ITYPE_t j) nogil:
        return self.distance_metric.rdist(&self.X[i, 0], &self.Y[j, 0], self.d)

    @final
    cdef DTYPE_t dist(self, ITYPE_t i, ITYPE_t j) nogil:
        return self.distance_metric.dist(&self.X[i, 0], &self.Y[j, 0], self.d)
