# Author: Carlos Perales <sir.perales@gmail.com>
# License: BSD 3 clause

from ..base import ClassifierMixin, BaseEstimator
from .kernel_ridge import _BaseKernelRidge


class KernelRidgeClassifier(_BaseKernelRidge, BaseEstimator, ClassifierMixin):
    """Kernel Ridge Classifier.

    Kernel Ridge Classifier (KRC) combines ridge
    regression (linear least squares with l2-norm regularization), special
    encoding for target and the kernel trick. It thus
    learns a linear function in the space induced by the respective kernel and
    the data. For non-linear kernels, this corresponds to a non-linear
    function in the original space.

    The form of the model learned by KRC is identical to support vector
    classification (SVC). However, different loss functions are used:
    KRR uses squared error loss while support vector classification uses
    epsilon-insensitive loss, both combined with l2 regularization.
    In contrast to SVC, fitting a multiclass KRC model can be done
    at once, without using one-against-one or one-against-rest techniques.

    Read more in the :ref:`User Guide <kernel_ridge>`.

    Parameters
    ----------
    alpha : {float, array-like}, shape = [n_targets]
        Small positive values of alpha improve the conditioning of the problem
        and reduce the variance of the estimates.  Alpha corresponds to
        ``(2*C)^-1`` in other linear models such as LogisticRegression or
        LinearSVC. If an array is passed, penalties are assumed to be specific
        to the targets. Hence they must correspond in number.

    kernel : string or callable, default="linear"
        Kernel mapping used internally. A callable should accept two arguments
        and the keyword arguments passed to this object as kernel_params, and
        should return a floating point number.

    gamma : float, default=None
        Gamma parameter for the RBF, laplacian, polynomial, exponential chi2
        and sigmoid kernels. Interpretation of the default value is left to
        the kernel; see the documentation for sklearn.metrics.pairwise.
        Ignored by other kernels.

    degree : float, default=3
        Degree of the polynomial kernel. Ignored by other kernels.

    coef0 : float, default=1
        Zero coefficient for polynomial and sigmoid kernels.
        Ignored by other kernels.

    kernel_params : mapping of string to any, optional
        Additional parameters (keyword arguments) for kernel function passed
        as callable object.

    Attributes
    ----------
    dual_coef_ : array, shape = [n_samples] or [n_samples, n_targets]
        Representation of weight vector(s) in kernel space

    X_fit_ : {array-like, sparse matrix}, shape = [n_samples, n_features]
        Training data, which is also required for prediction

    classes_ : array of shape = [n_classes] or a list of such arrays
        The classes labels (single output problem),
        or a list of arrays of class labels (multi-output problem).

    References
    ----------
    * Kevin P. Murphy
      "Machine Learning: A Probabilistic Perspective", The MIT Press
      chapter 14.4.3, pp. 492-493

    * An, Senjian, Wanquan Liu, and Svetha Venkatesh.
      "Face recognition using kernel ridge regression." Computer Vision
      and Pattern Recognition, 2007. CVPR'07. IEEE Conference on. IEEE,
      2007.

    See also
    --------
    sklearn.linear_model.Ridge:
        Linear ridge regression.
    sklearn.kernel_ridge.KernelRidge:
        Kernel Ridge implemented for regressions.

    """
    def __init__(self, alpha=0.1, kernel="rbf", gamma=None, degree=3, coef0=1,
                 kernel_params=None):
        self.alpha = alpha
        self.kernel = kernel
        self.gamma = gamma
        self.degree = degree
        self.coef0 = coef0
        self.kernel_params = kernel_params
