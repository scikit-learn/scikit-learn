import operator
import numpy as np


def l1_min_c(X, y, loss='l2', fit_intercept=True, intercept_scaling=1.0):
    """
    Return the maximum value for C that yields a model with coefficients
    and intercept set to zero for l1 penalized classifiers,
    such as LinearSVC with penalty='l1' and linear_model.LogisticRegression
    with penalty='l1'.

    This value is valid if class_weight parameter in fit() is not set.

    Parameters
    ----------
    X : array-like or sparse matrix, shape = [n_samples, n_features]
        Training vector, where n_samples in the number of samples and
        n_features is the number of features.

    y : array, shape = [n_samples]
        Target vector relative to X

    loss : {'l2', 'log'}, default to 'l2'
        Specifies the loss function.
        With 'l2' it is the l2 loss (a.k.a. squared hinge loss).
        With 'log' it is the loss of logistic regression models.

    fit_intercept : bool, default: True
        Specifies if the intercept should be fitted by the model.
        It must match the fit() method paramenter.

    intercept_scaling : float, default: 1
        when fit_intercept is True, instance vector x becomes
        [x, intercept_scaling],
        i.e. a "synthetic" feature with constant value equals to
        intercept_scaling is appended to the instance vector.
        It must match the fit() method parameter.

    Returns
    -------
    l1_min_c: float
        minimum value for C
    """
    import scipy.sparse as sp

    if loss not in ('l2', 'log'):
        raise ValueError('loss type not in ("l2", "log")')

    y = np.asarray(y)

    if sp.issparse(X):
        X = sp.csc_matrix(X)
        hstack = sp.hstack
        dot = operator.mul
    else:
        X = np.asarray(X)
        hstack = np.hstack
        dot = np.dot

    if fit_intercept:
        bias = intercept_scaling * np.ones((np.size(y), 1))
        X = hstack((X, bias))

    classes = np.unique(y)
    n_classes = np.size(classes)
    if n_classes <= 2:
        c = classes[0]
        y = y.reshape((1, -1))
        _y = np.empty(y.shape)
        _y[y == c] = 1
        _y[y != c] = -1
    else:
        _y = np.empty((n_classes, np.size(y)))
        for i, c in enumerate(classes):
            _y[i, y == c] = 1
            _y[i, y != c] = -1

    den = np.max(np.abs(dot(_y, X)))
    if den == 0.0:
        raise ValueError('Ill-posed l1_min_c calculation')
    if loss == 'l2':
        return 0.5 / den
    else:  # loss == 'log':
        return 2.0 / den
