
from ..base import ClassifierMixin, RegressorMixin
from .base import SparseBaseLibSVM


class SVC(SparseBaseLibSVM, ClassifierMixin):
    """SVC for sparse matrices (csr).

    See :class:`scikits.learn.svm.SVC` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).
    """

    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 cache_size=100.0, eps=1e-3, C=1.0, shrinking=True,
                 probability=False):

        SparseBaseLibSVM.__init__(self, 'c_svc', kernel, degree, gamma, coef0,
                         cache_size, eps, C, 0., 0.,
                         shrinking, probability)



class NuSVC (SparseBaseLibSVM, ClassifierMixin):
    """NuSVC for sparse matrices (csr).

    See :class:`scikits.learn.svm.NuSVC` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).
    """


    def __init__(self, nu=0.5, kernel='rbf', degree=3, gamma=0.0,
                 coef0=0.0, shrinking=True, probability=False,
                 eps=1e-3, cache_size=100.0):

        SparseBaseLibSVM.__init__(self, 'nu_svc', kernel, degree,
                         gamma, coef0, cache_size, eps, 0., nu, 0.,
                         shrinking, probability)




class SVR (SparseBaseLibSVM, RegressorMixin):
    """SVR for sparse matrices (csr)

    See :class:`scikits.learn.svm.SVR` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).
    """


    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 cache_size=100.0, eps=1e-3, C=1.0, nu=0.5, p=0.1,
                 shrinking=True, probability=False):

        SparseBaseLibSVM.__init__(self, 'epsilon_svr', kernel,
                         degree, gamma, coef0, cache_size, eps, C, nu,
                         p, shrinking, probability)





class NuSVR (SparseBaseLibSVM, RegressorMixin):
    """NuSVR for sparse matrices (csr)

    See :class:`scikits.learn.svm.NuSVC` for a complete list of parameters

    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).
    """

    def __init__(self, nu=0.5, C=1.0, kernel='rbf', degree=3,
                 gamma=0.0, coef0=0.0, shrinking=True,
                 probability=False, cache_size=100.0, eps=1e-3):

        SparseBaseLibSVM.__init__(self, 'epsilon_svr', kernel,
                         degree, gamma, coef0, cache_size, eps, C, nu,
                         0., shrinking, probability)



class OneClassSVM (SparseBaseLibSVM):
    """NuSVR for sparse matrices (csr)

    See :class:`scikits.learn.svm.NuSVC` for a complete list of parameters
    
    Notes
    -----
    For best results, this accepts a matrix in csr format
    (scipy.sparse.csr), but should be able to convert from any array-like
    object (including other sparse representations).
    """

    def __init__(self, kernel='rbf', degree=3, gamma=0.0, coef0=0.0,
                 cache_size=100.0, eps=1e-3, nu=0.5, shrinking=True,
                 probability=False):

        SparseBaseLibSVM.__init__(self, 'one_class', kernel, degree,
                         gamma, coef0, cache_size, eps, 0.0, nu, 0.0,
                         shrinking, probability)

    def fit(self, X, class_weight={}, sample_weight=[]):
        super(OneClassSVM, self).fit(
            X, [], class_weight=class_weight, ample_weight=sample_weight)
