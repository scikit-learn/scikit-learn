"""Partial dependence plots for tree ensembles. """

# Authors: Peter Prettenhofer
# License: BSD 3 clause

import warnings
from ..partial_dependence import partial_dependence as new_pd
from ..partial_dependence import plot_partial_dependence as new_ppd


def partial_dependence(gbrt, target_variables, grid=None, X=None,
                       percentiles=(0.05, 0.95), grid_resolution=100):
    """Partial dependence of ``target_variables``.

    Partial dependence plots show the dependence between the joint values
    of the ``target_variables`` and the function represented
    by the ``gbrt``.

    Read more in the :ref:`User Guide <partial_dependence>`.

    Parameters
    ----------
    gbrt : BaseGradientBoosting
        A fitted gradient boosting model.
    target_variables : array-like, dtype=int
        The target features for which the partial dependecy should be
        computed (size should be smaller than 3 for visual renderings).
    grid : array-like, shape=(n_points, len(target_variables))
        The grid of ``target_variables`` values for which the
        partial dependecy should be evaluated (either ``grid`` or ``X``
        must be specified).
    X : array-like, shape=(n_samples, n_features)
        The data on which ``gbrt`` was trained. It is used to generate
        a ``grid`` for the ``target_variables``. The ``grid`` comprises
        ``grid_resolution`` equally spaced points between the two
        ``percentiles``.
    percentiles : (low, high), default=(0.05, 0.95)
        The lower and upper percentile used create the extreme values
        for the ``grid``. Only if ``X`` is not None.
    grid_resolution : int, default=100
        The number of equally spaced points on the ``grid``.

    Returns
    -------
    pdp : array, shape=(n_classes, n_points)
        The partial dependence function evaluated on the ``grid``.
        For regression and binary classification ``n_classes==1``.
    axes : seq of ndarray or None
        The axes with which the grid has been created or None if
        the grid has been given.

    Examples
    --------
    >>> samples = [[0, 0, 2], [1, 0, 0]]
    >>> labels = [0, 1]
    >>> from sklearn.ensemble import GradientBoostingClassifier
    >>> gb = GradientBoostingClassifier(random_state=0).fit(samples, labels)
    >>> kwargs = dict(X=samples, percentiles=(0, 1), grid_resolution=2)
    >>> partial_dependence(gb, [0], **kwargs) # doctest: +SKIP
    (array([[-4.52...,  4.52...]]), [array([ 0.,  1.])])
    """
    warnings.warn("The function ensemble.partial_dependence has been moved to "
                  "partial_dependence in 0.20 and will be removed in 0.22.",
                  DeprecationWarning)
    return new_pd(est=gbrt,
                  target_variables=target_variables,
                  grid=grid,
                  X=X,
                  percentiles=percentiles,
                  grid_resolution=grid_resolution,
                  method='recursion')


def plot_partial_dependence(gbrt, X, features, feature_names=None,
                            label=None, n_cols=3, grid_resolution=100,
                            percentiles=(0.05, 0.95), n_jobs=None,
                            verbose=0, ax=None, line_kw=None,
                            contour_kw=None, **fig_kw):
    """Partial dependence plots for ``features``.

    The ``len(features)`` plots are arranged in a grid with ``n_cols``
    columns. Two-way partial dependence plots are plotted as contour
    plots.

    Read more in the :ref:`User Guide <partial_dependence>`.

    Parameters
    ----------
    gbrt : BaseGradientBoosting
        A fitted gradient boosting model.
    X : array-like, shape=(n_samples, n_features)
        The data on which ``gbrt`` was trained.
    features : seq of ints, strings, or tuples of ints or strings
        If seq[i] is an int or a tuple with one int value, a one-way
        PDP is created; if seq[i] is a tuple of two ints, a two-way
        PDP is created.
        If feature_names is specified and seq[i] is an int, seq[i]
        must be < len(feature_names).
        If seq[i] is a string, feature_names must be specified, and
        seq[i] must be in feature_names.
    feature_names : seq of str
        Name of each feature; feature_names[i] holds
        the name of the feature with index i.
    label : object
        The class label for which the PDPs should be computed.
        Only if gbrt is a multi-class model. Must be in ``gbrt.classes_``.
    n_cols : int
        The number of columns in the grid plot (default: 3).
    grid_resolution : int, default=100
        The number of equally spaced points on the axes.
    percentiles : (low, high), default=(0.05, 0.95)
        The lower and upper percentile used to create the extreme values
        for the PDP axes.
    n_jobs : int or None, optional (default=None)
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.
    verbose : int
        Verbose output during PD computations. Defaults to 0.
    ax : Matplotlib axis object, default None
        An axis object onto which the plots will be drawn.
    line_kw : dict
        Dict with keywords passed to the ``matplotlib.pyplot.plot`` call.
        For one-way partial dependence plots.
    contour_kw : dict
        Dict with keywords passed to the ``matplotlib.pyplot.plot`` call.
        For two-way partial dependence plots.
    **fig_kw : dict
        Dict with keywords passed to the figure() call.
        Note that all keywords not recognized above will be automatically
        included here.

    Returns
    -------
    fig : figure
        The Matplotlib Figure object.
    axs : seq of Axis objects
        A seq of Axis objects, one for each subplot.

    Examples
    --------
    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> X, y = make_friedman1()
    >>> clf = GradientBoostingRegressor(n_estimators=10).fit(X, y)
    >>> fig, axs = plot_partial_dependence(clf, X, [0, (0, 1)]) #doctest: +SKIP
    ...
    """
    warnings.warn("The function ensemble.plot_partial_dependence has been "
                  "moved to partial_dependence in 0.20 and will be removed "
                  "in 0.22.",
                  DeprecationWarning)
    return new_ppd(est=gbrt,
                   X=X,
                   features=features,
                   feature_names=feature_names,
                   label=label,
                   n_cols=n_cols,
                   grid_resolution=grid_resolution,
                   method='recursion',
                   percentiles=percentiles,
                   n_jobs=n_jobs,
                   verbose=verbose,
                   ax=ax,
                   line_kw=line_kw,
                   contour_kw=contour_kw,
                   **fig_kw)
