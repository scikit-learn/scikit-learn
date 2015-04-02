"""
Quadratic Discriminant Analysis
"""

# Author: Matthieu Perrot <matthieu.perrot@gmail.com>
#
# License: BSD 3 clause

import warnings

import numpy as np

from .discriminant_analysis import QuadraticDiscriminantAnalysis


class QDA(QuadraticDiscriminantAnalysis):
    """
    Quadratic Discriminant Analysis (QDA)

    A classifier with a quadratic decision boundary, generated
    by fitting class conditional densities to the data
    and using Bayes' rule.

    The model fits a Gaussian density to each class.

    Parameters
    ----------
    priors : array, optional, shape = [n_classes]
        Priors on classes

    reg_param : float, optional
        Regularizes the covariance estimate as
        ``(1-reg_param)*Sigma + reg_param*np.eye(n_features)``

    Attributes
    ----------
    covariances_ : list of array-like, shape = [n_features, n_features]
        Covariance matrices of each class.

    means_ : array-like, shape = [n_classes, n_features]
        Class means.

    priors_ : array-like, shape = [n_classes]
        Class priors (sum to 1).

    rotations_ : list of arrays
        For each class k an array of shape [n_features, n_k], with
        ``n_k = min(n_features, number of elements in class k)``
        It is the rotation of the Gaussian distribution, i.e. its
        principal axis.

    scalings_ : list of arrays
        For each class k an array of shape [n_k]. It contains the scaling
        of the Gaussian distributions along its principal axes, i.e. the
        variance in the rotated coordinate system.

    Examples
    --------
    >>> from sklearn.qda import QDA
    >>> import numpy as np
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> y = np.array([1, 1, 1, 2, 2, 2])
    >>> clf = QDA()
    >>> clf.fit(X, y)
    QDA(priors=None, reg_param=0.0)
    >>> print(clf.predict([[-0.8, -1]]))
    [1]

    See also
    --------
    sklearn.lda.LDA: Linear Discriminant Analysis
    """

    def __init__(self, priors=None, reg_param=0.):
        warnings.warn("qda.QDA has been moved to"
                      " discriminant_analysis.QuadraticDiscriminantAnalysis"
                      " in 0.17 and will be removed in 0.19.", DeprecationWarning)
        super(QDA, self).__init__(priors=priors, reg_param=reg_param)
