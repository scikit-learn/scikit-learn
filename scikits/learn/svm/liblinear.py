
from ..base import ClassifierMixin
from .base import BaseLibLinear


class LinearSVC(BaseLibLinear, ClassifierMixin):
    """
    Linear Support Vector Classification.

    Similar to SVC with parameter kernel='linear', but uses internally
    liblinear rather than libsvm, so it has more flexibility in the
    choice of penalties and loss functions and should be faster for
    huge datasets.

    Parameters
    ----------
    loss : string, 'l1' or 'l2' (default 'l2')
        Specifies the loss function. With 'l1' it is the standard SVM
        loss (a.k.a. hinge Loss) while with 'l2' it is the squared loss.
        (a.k.a. squared hinge Loss)

    penalty : string, 'l1' or 'l2' (default 'l2')
        Specifies the norm used in the penalization. The 'l2'
        penalty is the standard used in SVC. The 'l1' leads to coef_
        vectors that are sparse.

    dual : bool, (default True)
        Select the algorithm to either solve the dual or primal
        optimization problem.

    eps: float, optional
         precision for stopping criteria

    Attributes
    ----------
    `support_` : array-like, shape = [nSV, n_features]
        Support vectors.

    `dual_coef_` : array, shape = [n_class-1, nSV]
        Coefficients of the support vector in the decision function.

    `coef_` : array, shape = [n_class-1, n_features]
        Weights asigned to the features (coefficients in the primal
        problem). This is only available in the case of linear kernel.

    `intercept_` : array, shape = [n_class-1]
        Constants in decision function.

    Notes
    -----
    Some features of liblinear are still not wrapped, like the Cramer
    & Singer algorithm.

    """

    # all the implementation is provided by the mixins
    pass
