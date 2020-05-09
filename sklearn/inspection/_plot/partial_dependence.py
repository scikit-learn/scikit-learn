import numbers
from itertools import chain
from itertools import count
import warnings

import numpy as np
from scipy import sparse
from scipy.stats.mstats import mquantiles
from joblib import Parallel, delayed

from .. import partial_dependence
from ...base import is_regressor
from ...utils import check_array
from ...utils import check_matplotlib_support  # noqa
from ...utils import _safe_indexing
from ...utils.validation import _deprecate_positional_args


@_deprecate_positional_args
def plot_partial_dependence(estimator, X, features, *, feature_names=None,
                            target=None, response_method='auto', n_cols=3,
                            grid_resolution=100, percentiles=(0.05, 0.95),
                            method='auto', n_jobs=None, verbose=0, fig=None,
                            line_kw=None, contour_kw=None, ax=None):
    """Partial dependence plots.

    The ``len(features)`` plots are arranged in a grid with ``n_cols``
    columns. Two-way partial dependence plots are plotted as contour plots. The
    deciles of the feature values will be shown with tick marks on the x-axes
    for one-way plots, and on both axes for two-way plots.

    Read more in the :ref:`User Guide <partial_dependence>`.

    .. note::

        :func:`plot_partial_dependence` does not support using the same axes
        with multiple calls. To plot the the partial dependence for multiple
        estimators, please pass the axes created by the first call to the
        second call::

          >>> from sklearn.inspection import plot_partial_dependence
          >>> from sklearn.datasets import make_friedman1
          >>> from sklearn.linear_model import LinearRegression
          >>> X, y = make_friedman1()
          >>> est = LinearRegression().fit(X, y)
          >>> disp1 = plot_partial_dependence(est, X)  # doctest: +SKIP
          >>> disp2 = plot_partial_dependence(est, X,
          ...                                 ax=disp1.axes_)  # doctest: +SKIP

    .. warning::

        For :class:`~sklearn.ensemble.GradientBoostingClassifier` and
        :class:`~sklearn.ensemble.GradientBoostingRegressor`, the
        'recursion' method (used by default) will not account for the `init`
        predictor of the boosting process. In practice, this will produce
        the same values as 'brute' up to a constant offset in the target
        response, provided that `init` is a constant estimator (which is the
        default). However, if `init` is not a constant estimator, the
        partial dependence values are incorrect for 'recursion' because the
        offset will be sample-dependent. It is preferable to use the 'brute'
        method. Note that this only applies to
        :class:`~sklearn.ensemble.GradientBoostingClassifier` and
        :class:`~sklearn.ensemble.GradientBoostingRegressor`, not to
        :class:`~sklearn.ensemble.HistGradientBoostingClassifier` and
        :class:`~sklearn.ensemble.HistGradientBoostingRegressor`.

    Parameters
    ----------
    estimator : BaseEstimator
        A fitted estimator object implementing :term:`predict`,
        :term:`predict_proba`, or :term:`decision_function`.
        Multioutput-multiclass classifiers are not supported.

    X : {array-like or dataframe} of shape (n_samples, n_features)
        ``X`` is used to generate a grid of values for the target
        ``features`` (where the partial dependence will be evaluated), and
        also to generate values for the complement features when the
        `method` is 'brute'.

    features : list of {int, str, pair of int, pair of str}
        The target features for which to create the PDPs.
        If features[i] is an int or a string, a one-way PDP is created; if
        features[i] is a tuple, a two-way PDP is created. Each tuple must be
        of size 2.
        if any entry is a string, then it must be in ``feature_names``.

    feature_names : array-like of shape (n_features,), dtype=str, default=None
        Name of each feature; feature_names[i] holds the name of the feature
        with index i.
        By default, the name of the feature corresponds to their numerical
        index for NumPy array and their column name for pandas dataframe.

    target : int, optional (default=None)
        - In a multiclass setting, specifies the class for which the PDPs
          should be computed. Note that for binary classification, the
          positive class (index 1) is always used.
        - In a multioutput setting, specifies the task for which the PDPs
          should be computed.

        Ignored in binary classification or classical regression settings.

    response_method : 'auto', 'predict_proba' or 'decision_function', \
            optional (default='auto')
        Specifies whether to use :term:`predict_proba` or
        :term:`decision_function` as the target response. For regressors
        this parameter is ignored and the response is always the output of
        :term:`predict`. By default, :term:`predict_proba` is tried first
        and we revert to :term:`decision_function` if it doesn't exist. If
        ``method`` is 'recursion', the response is always the output of
        :term:`decision_function`.

    n_cols : int, optional (default=3)
        The maximum number of columns in the grid plot. Only active when `ax`
        is a single axis or `None`.

    grid_resolution : int, optional (default=100)
        The number of equally spaced points on the axes of the plots, for each
        target feature.

    percentiles : tuple of float, optional (default=(0.05, 0.95))
        The lower and upper percentile used to create the extreme values
        for the PDP axes. Must be in [0, 1].

    method : str, optional (default='auto')
        The method used to calculate the averaged predictions:

        - 'recursion' is only supported for some tree-based estimators (namely
          :class:`~sklearn.ensemble.GradientBoostingClassifier`,
          :class:`~sklearn.ensemble.GradientBoostingRegressor`,
          :class:`~sklearn.ensemble.HistGradientBoostingClassifier`,
          :class:`~sklearn.ensemble.HistGradientBoostingRegressor`,
          :class:`~sklearn.tree.DecisionTreeRegressor`,
          :class:`~sklearn.ensemble.RandomForestRegressor`
          but is more efficient in terms of speed.
          With this method, the target response of a
          classifier is always the decision function, not the predicted
          probabilities.

        - 'brute' is supported for any estimator, but is more
          computationally intensive.

        - 'auto': the 'recursion' is used for estimators that support it,
          and 'brute' is used otherwise.

        Please see :ref:`this note <pdp_method_differences>` for
        differences between the 'brute' and 'recursion' method.

    n_jobs : int, optional (default=None)
        The number of CPUs to use to compute the partial dependences.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.

    verbose : int, optional (default=0)
        Verbose output during PD computations.

    fig : Matplotlib figure object, optional (default=None)
        A figure object onto which the plots will be drawn, after the figure
        has been cleared. By default, a new one is created.

        .. deprecated:: 0.22
           ``fig`` will be removed in 0.24.

    line_kw : dict, optional
        Dict with keywords passed to the ``matplotlib.pyplot.plot`` call.
        For one-way partial dependence plots.

    contour_kw : dict, optional
        Dict with keywords passed to the ``matplotlib.pyplot.contourf`` call.
        For two-way partial dependence plots.

    ax : Matplotlib axes or array-like of Matplotlib axes, default=None
        - If a single axis is passed in, it is treated as a bounding axes
            and a grid of partial dependence plots will be drawn within
            these bounds. The `n_cols` parameter controls the number of
            columns in the grid.
        - If an array-like of axes are passed in, the partial dependence
            plots will be drawn directly into these axes.
        - If `None`, a figure and a bounding axes is created and treated
            as the single axes case.

        .. versionadded:: 0.22

    Returns
    -------
    display: :class:`~sklearn.inspection.PartialDependenceDisplay`

    Examples
    --------
    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> X, y = make_friedman1()
    >>> clf = GradientBoostingRegressor(n_estimators=10).fit(X, y)
    >>> plot_partial_dependence(clf, X, [0, (0, 1)]) #doctest: +SKIP

    See also
    --------
    sklearn.inspection.partial_dependence: Return raw partial
      dependence values
    """
    check_matplotlib_support('plot_partial_dependence')  # noqa
    import matplotlib.pyplot as plt  # noqa
    from matplotlib import transforms  # noqa
    from matplotlib.ticker import MaxNLocator  # noqa
    from matplotlib.ticker import ScalarFormatter  # noqa

    # set target_idx for multi-class estimators
    if hasattr(estimator, 'classes_') and np.size(estimator.classes_) > 2:
        if target is None:
            raise ValueError('target must be specified for multi-class')
        target_idx = np.searchsorted(estimator.classes_, target)
        if (not (0 <= target_idx < len(estimator.classes_)) or
                estimator.classes_[target_idx] != target):
            raise ValueError('target not in est.classes_, got {}'.format(
                target))
    else:
        # regression and binary classification
        target_idx = 0

    # Use check_array only on lists and other non-array-likes / sparse. Do not
    # convert DataFrame into a NumPy array.
    if not(hasattr(X, '__array__') or sparse.issparse(X)):
        X = check_array(X, force_all_finite='allow-nan', dtype=np.object)
    n_features = X.shape[1]

    # convert feature_names to list
    if feature_names is None:
        if hasattr(X, "loc"):
            # get the column names for a pandas dataframe
            feature_names = X.columns.tolist()
        else:
            # define a list of numbered indices for a numpy array
            feature_names = [str(i) for i in range(n_features)]
    elif hasattr(feature_names, "tolist"):
        # convert numpy array or pandas index to a list
        feature_names = feature_names.tolist()
    if len(set(feature_names)) != len(feature_names):
        raise ValueError('feature_names should not contain duplicates.')

    def convert_feature(fx):
        if isinstance(fx, str):
            try:
                fx = feature_names.index(fx)
            except ValueError:
                raise ValueError('Feature %s not in feature_names' % fx)
        return int(fx)

    # convert features into a seq of int tuples
    tmp_features = []
    for fxs in features:
        if isinstance(fxs, (numbers.Integral, str)):
            fxs = (fxs,)
        try:
            fxs = tuple(convert_feature(fx) for fx in fxs)
        except TypeError:
            raise ValueError('Each entry in features must be either an int, '
                             'a string, or an iterable of size at most 2.')
        if not 1 <= np.size(fxs) <= 2:
            raise ValueError('Each entry in features must be either an int, '
                             'a string, or an iterable of size at most 2.')

        tmp_features.append(fxs)

    features = tmp_features

    # Early exit if the axes does not have the correct number of axes
    if ax is not None and not isinstance(ax, plt.Axes):
        axes = np.asarray(ax, dtype=object)
        if axes.size != len(features):
            raise ValueError("Expected ax to have {} axes, got {}".format(
                             len(features), axes.size))

    for i in chain.from_iterable(features):
        if i >= len(feature_names):
            raise ValueError('All entries of features must be less than '
                             'len(feature_names) = {0}, got {1}.'
                             .format(len(feature_names), i))

    # compute averaged predictions
    pd_results = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(partial_dependence)(estimator, X, fxs,
                                    response_method=response_method,
                                    method=method,
                                    grid_resolution=grid_resolution,
                                    percentiles=percentiles)
        for fxs in features)

    # For multioutput regression, we can only check the validity of target
    # now that we have the predictions.
    # Also note: as multiclass-multioutput classifiers are not supported,
    # multiclass and multioutput scenario are mutually exclusive. So there is
    # no risk of overwriting target_idx here.
    avg_preds, _ = pd_results[0]  # checking the first result is enough
    if is_regressor(estimator) and avg_preds.shape[0] > 1:
        if target is None:
            raise ValueError(
                'target must be specified for multi-output regressors')
        if not 0 <= target <= avg_preds.shape[0]:
            raise ValueError(
                'target must be in [0, n_tasks], got {}.'.format(target))
        target_idx = target

    # get global min and max average predictions of PD grouped by plot type
    pdp_lim = {}
    for avg_preds, values in pd_results:
        min_pd = avg_preds[target_idx].min()
        max_pd = avg_preds[target_idx].max()
        n_fx = len(values)
        old_min_pd, old_max_pd = pdp_lim.get(n_fx, (min_pd, max_pd))
        min_pd = min(min_pd, old_min_pd)
        max_pd = max(max_pd, old_max_pd)
        pdp_lim[n_fx] = (min_pd, max_pd)

    deciles = {}
    for fx in chain.from_iterable(features):
        if fx not in deciles:
            X_col = _safe_indexing(X, fx, axis=1)
            deciles[fx] = mquantiles(X_col, prob=np.arange(0.1, 1.0, 0.1))

    if fig is not None:
        warnings.warn("The fig parameter is deprecated in version "
                      "0.22 and will be removed in version 0.24",
                      FutureWarning)
        fig.clear()
        ax = fig.gca()

    display = PartialDependenceDisplay(pd_results=pd_results,
                                       features=features,
                                       feature_names=feature_names,
                                       target_idx=target_idx,
                                       pdp_lim=pdp_lim,
                                       deciles=deciles)
    return display.plot(ax=ax, n_cols=n_cols, line_kw=line_kw,
                        contour_kw=contour_kw)


class PartialDependenceDisplay:
    """Partial Dependence Plot (PDP) visualization.

    It is recommended to use
    :func:`~sklearn.inspection.plot_partial_dependence` to create a
    :class:`~sklearn.inspection.PartialDependenceDisplay`. All parameters are
    stored as attributes.

    Read more in
    :ref:`sphx_glr_auto_examples_miscellaneous_plot_partial_dependence_visualization_api.py`
    and the :ref:`User Guide <visualizations>`.

        .. versionadded:: 0.22

    Parameters
    ----------
    pd_results : list of (ndarray, ndarray)
        Results of :func:`~sklearn.inspection.partial_dependence` for
        ``features``. Each tuple corresponds to a (averaged_predictions, grid).

    features : list of (int,) or list of (int, int)
        Indices of features for a given plot. A tuple of one integer will plot
        a partial dependence curve of one feature. A tuple of two integers will
        plot a two-way partial dependence curve as a contour plot.

    feature_names : list of str
        Feature names corresponding to the indices in ``features``.

    target_idx : int

        - In a multiclass setting, specifies the class for which the PDPs
          should be computed. Note that for binary classification, the
          positive class (index 1) is always used.
        - In a multioutput setting, specifies the task for which the PDPs
          should be computed.

        Ignored in binary classification or classical regression settings.

    pdp_lim : dict
        Global min and max average predictions, such that all plots will have
        the same scale and y limits. `pdp_lim[1]` is the global min and max for
        single partial dependence curves. `pdp_lim[2]` is the global min and
        max for two-way partial dependence curves.

    deciles : dict
        Deciles for feature indices in ``features``.

    Attributes
    ----------
    bounding_ax_ : matplotlib Axes or None
        If `ax` is an axes or None, the `bounding_ax_` is the axes where the
        grid of partial dependence plots are drawn. If `ax` is a list of axes
        or a numpy array of axes, `bounding_ax_` is None.

    axes_ : ndarray of matplotlib Axes
        If `ax` is an axes or None, `axes_[i, j]` is the axes on the i-th row
        and j-th column. If `ax` is a list of axes, `axes_[i]` is the i-th item
        in `ax`. Elements that are None correspond to a nonexisting axes in
        that position.

    lines_ : ndarray of matplotlib Artists
        If `ax` is an axes or None, `lines_[i, j]` is the partial dependence
        curve on the i-th row and j-th column. If `ax` is a list of axes,
        `lines_[i]` is the partial dependence curve corresponding to the i-th
        item in `ax`. Elements that are None correspond to a nonexisting axes
        or an axes that does not include a line plot.

    deciles_vlines_ : ndarray of matplotlib LineCollection
        If `ax` is an axes or None, `vlines_[i, j]` is the line collection
        representing the x axis deciles of the i-th row and j-th column. If
        `ax` is a list of axes, `vlines_[i]` corresponds to the i-th item in
        `ax`. Elements that are None correspond to a nonexisting axes or an
        axes that does not include a PDP plot.
        .. versionadded:: 0.23
    deciles_hlines_ : ndarray of matplotlib LineCollection
        If `ax` is an axes or None, `vlines_[i, j]` is the line collection
        representing the y axis deciles of the i-th row and j-th column. If
        `ax` is a list of axes, `vlines_[i]` corresponds to the i-th item in
        `ax`. Elements that are None correspond to a nonexisting axes or an
        axes that does not include a 2-way plot.
        .. versionadded:: 0.23

    contours_ : ndarray of matplotlib Artists
        If `ax` is an axes or None, `contours_[i, j]` is the partial dependence
        plot on the i-th row and j-th column. If `ax` is a list of axes,
        `contours_[i]` is the partial dependence plot corresponding to the i-th
        item in `ax`. Elements that are None correspond to a nonexisting axes
        or an axes that does not include a contour plot.

    figure_ : matplotlib Figure
        Figure containing partial dependence plots.

    """
    @_deprecate_positional_args
    def __init__(self, pd_results, *, features, feature_names, target_idx,
                 pdp_lim, deciles):
        self.pd_results = pd_results
        self.features = features
        self.feature_names = feature_names
        self.target_idx = target_idx
        self.pdp_lim = pdp_lim
        self.deciles = deciles

    def plot(self, ax=None, n_cols=3, line_kw=None, contour_kw=None):
        """Plot partial dependence plots.

        Parameters
        ----------
        ax : Matplotlib axes or array-like of Matplotlib axes, default=None
            - If a single axis is passed in, it is treated as a bounding axes
                and a grid of partial dependence plots will be drawn within
                these bounds. The `n_cols` parameter controls the number of
                columns in the grid.
            - If an array-like of axes are passed in, the partial dependence
                plots will be drawn directly into these axes.
            - If `None`, a figure and a bounding axes is created and treated
                as the single axes case.

        n_cols : int, default=3
            The maximum number of columns in the grid plot. Only active when
            `ax` is a single axes or `None`.

        line_kw : dict, default=None
            Dict with keywords passed to the `matplotlib.pyplot.plot` call.
            For one-way partial dependence plots.

        contour_kw : dict, default=None
            Dict with keywords passed to the `matplotlib.pyplot.contourf`
            call for two-way partial dependence plots.

        Returns
        -------
        display: :class:`~sklearn.inspection.PartialDependenceDisplay`
        """

        check_matplotlib_support("plot_partial_dependence")
        import matplotlib.pyplot as plt  # noqa
        from matplotlib import transforms  # noqa
        from matplotlib.ticker import MaxNLocator  # noqa
        from matplotlib.ticker import ScalarFormatter  # noqa
        from matplotlib.gridspec import GridSpecFromSubplotSpec  # noqa

        if line_kw is None:
            line_kw = {}
        if contour_kw is None:
            contour_kw = {}

        if ax is None:
            _, ax = plt.subplots()

        default_contour_kws = {"alpha": 0.75}
        contour_kw = {**default_contour_kws, **contour_kw}

        n_features = len(self.features)

        if isinstance(ax, plt.Axes):
            # If ax was set off, it has most likely been set to off
            # by a previous call to plot.
            if not ax.axison:
                raise ValueError("The ax was already used in another plot "
                                 "function, please set ax=display.axes_ "
                                 "instead")

            ax.set_axis_off()
            self.bounding_ax_ = ax
            self.figure_ = ax.figure

            n_cols = min(n_cols, n_features)
            n_rows = int(np.ceil(n_features / float(n_cols)))

            self.axes_ = np.empty((n_rows, n_cols), dtype=np.object)

            axes_ravel = self.axes_.ravel()

            gs = GridSpecFromSubplotSpec(n_rows, n_cols,
                                         subplot_spec=ax.get_subplotspec())
            for i, spec in zip(range(n_features), gs):
                axes_ravel[i] = self.figure_.add_subplot(spec)

        else:  # array-like
            ax = np.asarray(ax, dtype=object)
            if ax.size != n_features:
                raise ValueError("Expected ax to have {} axes, got {}"
                                 .format(n_features, ax.size))

            if ax.ndim == 2:
                n_cols = ax.shape[1]
            else:
                n_cols = None

            self.bounding_ax_ = None
            self.figure_ = ax.ravel()[0].figure
            self.axes_ = ax

        # create contour levels for two-way plots
        if 2 in self.pdp_lim:
            Z_level = np.linspace(*self.pdp_lim[2], num=8)

        self.lines_ = np.empty_like(self.axes_, dtype=np.object)
        self.contours_ = np.empty_like(self.axes_, dtype=np.object)
        self.deciles_vlines_ = np.empty_like(self.axes_, dtype=np.object)
        self.deciles_hlines_ = np.empty_like(self.axes_, dtype=np.object)
        # Create 1d views of these 2d arrays for easy indexing
        lines_ravel = self.lines_.ravel(order='C')
        contours_ravel = self.contours_.ravel(order='C')
        vlines_ravel = self.deciles_vlines_.ravel(order='C')
        hlines_ravel = self.deciles_hlines_.ravel(order='C')

        for i, axi, fx, (avg_preds, values) in zip(count(),
                                                   self.axes_.ravel(),
                                                   self.features,
                                                   self.pd_results):
            if len(values) == 1:
                lines_ravel[i] = axi.plot(values[0],
                                          avg_preds[self.target_idx].ravel(),
                                          **line_kw)[0]
            else:
                # contour plot
                XX, YY = np.meshgrid(values[0], values[1])
                Z = avg_preds[self.target_idx].T
                CS = axi.contour(XX, YY, Z, levels=Z_level, linewidths=0.5,
                                 colors='k')
                contours_ravel[i] = axi.contourf(XX, YY, Z, levels=Z_level,
                                                 vmax=Z_level[-1],
                                                 vmin=Z_level[0],
                                                 **contour_kw)
                axi.clabel(CS, fmt='%2.2f', colors='k', fontsize=10,
                           inline=True)

            trans = transforms.blended_transform_factory(axi.transData,
                                                         axi.transAxes)
            ylim = axi.get_ylim()
            vlines_ravel[i] = axi.vlines(self.deciles[fx[0]], 0, 0.05,
                                         transform=trans, color='k')
            axi.set_ylim(ylim)

            # Set xlabel if it is not already set
            if not axi.get_xlabel():
                axi.set_xlabel(self.feature_names[fx[0]])

            if len(values) == 1:
                if n_cols is None or i % n_cols == 0:
                    axi.set_ylabel('Partial dependence')
                else:
                    axi.set_yticklabels([])
                axi.set_ylim(self.pdp_lim[1])
            else:
                # contour plot
                trans = transforms.blended_transform_factory(axi.transAxes,
                                                             axi.transData)
                xlim = axi.get_xlim()
                hlines_ravel[i] = axi.hlines(self.deciles[fx[1]], 0, 0.05,
                                             transform=trans, color='k')
                # hline erases xlim
                axi.set_ylabel(self.feature_names[fx[1]])
                axi.set_xlim(xlim)
        return self
