from .. import LinearSVC
from ..base import BaseSVC
from ...base import RegressorMixin
from .base import SparseBaseLibSVM
from ...utils import deprecated


@deprecated("""to be removed in v0.14;
use sklearn.svm.SVC instead""")
class SVC(SparseBaseLibSVM, BaseSVC):
    """SVC for sparse matrices (csr).

    See :class:`sklearn.svm.SVC` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.svm.sparse import SVC
    >>> clf = SVC()
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
            gamma=0.0, kernel='rbf', probability=False, shrinking=True,
            tol=0.001, verbose=False)
    >>> print(clf.predict([[-0.8, -1]]))
    [ 1.]
    """

    def __init__(self, C=1.0, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=200, class_weight=None,
                 verbose=False):

        super(SVC, self).__init__('c_svc', kernel, degree, gamma, coef0, tol,
                                  C, 0., 0., shrinking, probability,
                                  cache_size, class_weight, verbose)


@deprecated("""to be removed in v0.14;
use sklearn.svm.NuSVC instead""")
class NuSVC(SparseBaseLibSVM, BaseSVC):
    """NuSVC for sparse matrices (csr).

    See :class:`sklearn.svm.NuSVC` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).

    Examples
    --------
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [1, 1], [2, 1]])
    >>> y = np.array([1, 1, 2, 2])
    >>> from sklearn.svm.sparse import NuSVC
    >>> clf = NuSVC()
    >>> clf.fit(X, y) #doctest: +NORMALIZE_WHITESPACE
    NuSVC(cache_size=200, class_weight=None, coef0=0.0, degree=3, gamma=0.0,
            kernel='rbf', nu=0.5, probability=False, shrinking=True, tol=0.001,
            verbose=False)
    >>> print(clf.predict([[-0.8, -1]]))
    [ 1.]
    """

    def __init__(self, nu=0.5, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 tol=1e-3, cache_size=200, class_weight=None,
                 verbose=False):

        super(NuSVC, self).__init__('nu_svc', kernel, degree, gamma, coef0,
                                    tol, 0., nu, 0., shrinking, probability,
                                    cache_size, class_weight, verbose)


@deprecated("""to be removed in v0.14;
use sklearn.svm.SVR instead""")
class SVR(SparseBaseLibSVM, RegressorMixin):
    """SVR for sparse matrices (csr)

    See :class:`sklearn.svm.SVR` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).

    Examples
    --------
    >>> from sklearn.svm.sparse import SVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = SVR(C=1.0, epsilon=0.2)
    >>> clf.fit(X, y)
    SVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.2, gamma=0.0,
      kernel='rbf', probability=False, shrinking=True, tol=0.001,
      verbose=False)
    """

    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0, tol=1e-3,
            C=1.0, epsilon=0.1, shrinking=True, probability=False,
            cache_size=200, verbose=False):

        super(SVR, self).__init__('epsilon_svr', kernel, degree, gamma, coef0,
                                  tol, C, 0., epsilon, shrinking, probability,
                                  cache_size, None, verbose)


@deprecated("""to be removed in v0.14;
use sklearn.svm.NuSVR instead""")
class NuSVR(SparseBaseLibSVM, RegressorMixin):
    """NuSVR for sparse matrices (csr)

    See :class:`sklearn.svm.NuSVC` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).

    Examples
    --------
    >>> from sklearn.svm.sparse import NuSVR
    >>> import numpy as np
    >>> n_samples, n_features = 10, 5
    >>> np.random.seed(0)
    >>> y = np.random.randn(n_samples)
    >>> X = np.random.randn(n_samples, n_features)
    >>> clf = NuSVR(nu=0.1, C=1.0)
    >>> clf.fit(X, y)
    NuSVR(C=1.0, cache_size=200, coef0=0.0, degree=3, epsilon=0.1, gamma=0.0,
       kernel='rbf', nu=0.1, probability=False, shrinking=True, tol=0.001,
       verbose=False)
    """

    def __init__(self, nu=0.5, C=1.0, kernel='rbf', degree=3, gamma=0.0,
            coef0=0.0, shrinking=True, epsilon=0.1, probability=False,
            tol=1e-3, cache_size=200, verbose=False):

        super(NuSVR, self).__init__('nu_svr', kernel, degree, gamma, coef0,
                tol, C, nu, epsilon, shrinking, probability, cache_size,
                None, verbose)


@deprecated("""to be removed in v0.14;
use sklearn.svm.OneClassSVM instead""")
class OneClassSVM(SparseBaseLibSVM):
    """OneClassSVM for sparse matrices (csr)

    See :class:`sklearn.svm.OneClassSVM` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).
    """

    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0, tol=1e-3,
            nu=0.5, shrinking=True, probability=False, cache_size=200,
            verbose=False):

        super(OneClassSVM, self).__init__('one_class', kernel, degree, gamma,
                coef0, tol, 0.0, nu, 0.0, shrinking, probability, cache_size,
                None, verbose)

    def fit(self, X, sample_weight=None):
        super(OneClassSVM, self).fit(
            X, [], sample_weight=sample_weight)


@deprecated("""to be removed in v0.14;
use sklearn.svm.LinearSVC instead""")
class LinearSVC(LinearSVC):
    pass
