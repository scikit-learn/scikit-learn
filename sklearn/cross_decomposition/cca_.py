from .pls_ import _PLS

__all__ = ['CCA']


class CCA(_PLS):
    """CCA Canonical Correlation Analysis.

    CCA inherits from PLS with mode="B" and deflation_mode="canonical".

    Read more in the :ref:`User Guide <cross_decomposition>`.

    Parameters
    ----------
    n_components : int, (default 2).
        number of components to keep.

    scale : boolean, (default True)
        whether to scale the data?

    max_iter : an integer, (default 500)
        the maximum number of iterations of the NIPALS inner loop

    tol : non-negative real, default 1e-06.
        the tolerance used in the iterative algorithm

    copy : boolean
        Whether the deflation be done on a copy. Let the default value
        to True unless you don't care about side effects

    Attributes
    ----------
    x_weights_ : array, [p, n_components]
        X block weights vectors.

    y_weights_ : array, [q, n_components]
        Y block weights vectors.

    x_loadings_ : array, [p, n_components]
        X block loadings vectors.

    y_loadings_ : array, [q, n_components]
        Y block loadings vectors.

    x_scores_ : array, [n_samples, n_components]
        X scores.

    y_scores_ : array, [n_samples, n_components]
        Y scores.

    x_rotations_ : array, [p, n_components]
        X block to latents rotations.

    y_rotations_ : array, [q, n_components]
        Y block to latents rotations.

    n_iter_ : array-like
        Number of iterations of the NIPALS inner loop for each
        component.

    Notes
    -----
    For each component k, find the weights u, v that maximizes
    max corr(Xk u, Yk v), such that ``|u| = |v| = 1``

    Note that it maximizes only the correlations between the scores.

    The residual matrix of X (Xk+1) block is obtained by the deflation on the
    current X score: x_score.

    The residual matrix of Y (Yk+1) block is obtained by deflation on the
    current Y score.

    Examples
    --------
    >>> from sklearn.cross_decomposition import CCA
    >>> X = [[0., 0., 1.], [1.,0.,0.], [2.,2.,2.], [3.,5.,4.]]
    >>> Y = [[0.1, -0.2], [0.9, 1.1], [6.2, 5.9], [11.9, 12.3]]
    >>> cca = CCA(n_components=1)
    >>> cca.fit(X, Y)
    ... # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
    CCA(copy=True, max_iter=500, n_components=1, scale=True, tol=1e-06)
    >>> X_c, Y_c = cca.transform(X, Y)

    References
    ----------

    Jacob A. Wegelin. A survey of Partial Least Squares (PLS) methods, with
    emphasis on the two-block case. Technical Report 371, Department of
    Statistics, University of Washington, Seattle, 2000.

    In french but still a reference:
    Tenenhaus, M. (1998). La regression PLS: theorie et pratique. Paris:
    Editions Technic.

    See also
    --------
    PLSCanonical
    PLSSVD
    """

    def __init__(self, n_components=2, scale=True,
                 max_iter=500, tol=1e-06, copy=True):
        super(CCA, self).__init__(n_components=n_components, scale=scale,
                                  deflation_mode="canonical", mode="B",
                                  norm_y_weights=True, algorithm="nipals",
                                  max_iter=max_iter, tol=tol, copy=copy)
