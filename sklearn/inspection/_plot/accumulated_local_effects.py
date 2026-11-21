# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import numbers
from math import ceil

import numpy as np

from sklearn.base import is_regressor
from sklearn.inspection._accumulated_local_effects import accumulated_local_effects
from sklearn.inspection._pd_utils import _check_feature_names, _get_feature_index
from sklearn.utils import check_array
from sklearn.utils._optional_dependencies import check_matplotlib_support
from sklearn.utils.parallel import Parallel, delayed


class AccumulatedLocalEffectsDisplay:
    """Accumulated Local Effects (ALE) visualization.

    It is recommended to use
    :func:`~sklearn.inspection.AccumulatedLocalEffectsDisplay.from_estimator` to
    create an :class:`~sklearn.inspection.AccumulatedLocalEffectsDisplay`. All
    parameters are stored as attributes.

    For general information regarding `scikit-learn` visualization tools, see
    the :ref:`Visualization Guide <visualizations>`.
    For guidance on interpreting these plots, refer to the
    :ref:`Inspection Guide <accumulated_local_effects>`.

    .. versionadded:: 1.10

    Parameters
    ----------
    ale_results : list of Bunch
        Results of :func:`~sklearn.inspection.accumulated_local_effects` for
        each entry of `features`.

    features : list of (int,) or list of (int, int)
        Indices of features for a given plot. A tuple of one integer plots the
        accumulated local effects of one feature, as a line (continuous) or bars
        (categorical). A tuple of two integers plots the second-order effects of
        two features as a heatmap.

    feature_names : list of str
        Feature names corresponding to the indices in `features`.

    target_idx : int
        The index of the class (multi-class) or task (multi-output regression)
        to plot. Ignored for binary classification and single-output
        regression.

    Attributes
    ----------
    axes_ : ndarray of matplotlib Axes
        The axes of the plots.

    lines_ : ndarray of matplotlib Line2D
        The lines of the continuous features. `None` for categorical features.

    bars_ : ndarray of matplotlib BarContainer
        The bars of the categorical features. `None` for continuous features.

    heatmaps_ : ndarray of matplotlib QuadMesh
        The heatmaps of the pairs of features. `None` for single features.

    figure_ : matplotlib Figure
        Figure containing the plots.

    See Also
    --------
    accumulated_local_effects : Compute accumulated local effects.
    PartialDependenceDisplay : Partial dependence visualization.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> from sklearn.datasets import make_friedman1
    >>> from sklearn.ensemble import GradientBoostingRegressor
    >>> from sklearn.inspection import AccumulatedLocalEffectsDisplay
    >>> X, y = make_friedman1(random_state=0)
    >>> reg = GradientBoostingRegressor(n_estimators=10).fit(X, y)
    >>> AccumulatedLocalEffectsDisplay.from_estimator(reg, X, [0, 3, (0, 1)])
    <...>
    >>> plt.show()
    """

    def __init__(self, ale_results, *, features, feature_names, target_idx):
        self.ale_results = ale_results
        self.features = features
        self.feature_names = feature_names
        self.target_idx = target_idx

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        features,
        *,
        categorical_features="from_dtype",
        feature_names=None,
        target=None,
        response_method="auto",
        grid_resolution=100,
        n_cols=3,
        n_jobs=None,
        ax=None,
        line_kw=None,
        bar_kw=None,
        heatmap_kw=None,
    ):
        """Plot accumulated local effects from an estimator.

        Read more in the :ref:`User Guide <accumulated_local_effects>`.

        .. versionadded:: 1.10

        Parameters
        ----------
        estimator : BaseEstimator
            A fitted estimator object implementing :term:`predict`,
            :term:`predict_proba`, or :term:`decision_function`.
            Multioutput-multiclass classifiers are not supported.

        X : {array-like, dataframe} of shape (n_samples, n_features)
            Data used to compute the accumulated local effects. Samples with a
            missing value in a plotted feature are ignored for that plot.

        features : list of {int, str, pair of int, pair of str}
            The target features for which to create the plots. A single feature
            gives a line (continuous) or bar (categorical) plot. A pair of
            continuous features, e.g. `(0, 1)`, gives a heatmap of their
            second-order effects, where cells without samples are left blank.

        categorical_features : array-like of {bool, int, str} of shape \
                (n_features,) or shape (n_categorical_features,), "from_dtype" \
                or None, default="from_dtype"
            Indicates the categorical or discrete features. See
            :func:`~sklearn.inspection.accumulated_local_effects`.

        feature_names : array-like of shape (n_features,), dtype=str, default=None
            Name of each feature; `feature_names[i]` holds the name of the
            feature with index `i`.

        target : int, default=None
            - In a multiclass setting, specifies the class to plot.
            - In a multioutput regression setting, specifies the task to plot.

            Ignored in binary classification or classical regression settings.

        response_method : {'auto', 'predict_proba', 'decision_function'}, \
                default='auto'
            Specifies whether to use :term:`predict_proba` or
            :term:`decision_function` as the target response.

        grid_resolution : int, default=100
            Maximum number of points at which the accumulated local effects of a
            continuous feature are computed. The range of the feature is split
            into `grid_resolution - 1` intervals with roughly the same number of
            samples.

        n_cols : int, default=3
            The maximum number of columns in the grid plot. Only active when
            `ax` is `None`.

        n_jobs : int, default=None
            The number of CPUs to use to compute the accumulated local effects.
            ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
            ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
            for more details.

        ax : Matplotlib axes or array-like of Matplotlib axes, default=None
            Axes to plot on, one per feature. If `None`, a figure is created.

        line_kw : dict, default=None
            Keyword arguments passed to ``matplotlib.pyplot.plot``.

        bar_kw : dict, default=None
            Keyword arguments passed to ``matplotlib.pyplot.bar``.

        heatmap_kw : dict, default=None
            Keyword arguments passed to ``matplotlib.pyplot.pcolormesh`` for the
            pairs of features.

        Returns
        -------
        display : :class:`~sklearn.inspection.AccumulatedLocalEffectsDisplay`
            Object that stores the computed values and the plots.
        """
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        if hasattr(estimator, "classes_") and np.size(estimator.classes_) > 2:
            if target is None:
                raise ValueError("target must be specified for multi-class")
            target_idx = np.searchsorted(estimator.classes_, target)
            if (
                not (0 <= target_idx < len(estimator.classes_))
                or estimator.classes_[target_idx] != target
            ):
                raise ValueError(f"target not in est.classes_, got {target}")
        else:
            target_idx = 0

        if not hasattr(X, "__array__"):
            X = check_array(X, ensure_all_finite="allow-nan", dtype=object)
        feature_names = _check_feature_names(X, feature_names)

        tmp_features = []
        for fxs in features:
            if isinstance(fxs, (numbers.Integral, str)):
                fxs = (fxs,)
            fxs = tuple(
                _get_feature_index(fx, feature_names=feature_names) for fx in fxs
            )
            if len(fxs) not in (1, 2):
                raise ValueError(
                    "Each entry in features must be an int, a string, or a pair of "
                    "them."
                )
            tmp_features.append(fxs)
        features = tmp_features

        ale_results = Parallel(n_jobs=n_jobs)(
            delayed(accumulated_local_effects)(
                estimator,
                X,
                fxs,
                categorical_features=categorical_features,
                feature_names=feature_names,
                response_method=response_method,
                grid_resolution=grid_resolution,
            )
            for fxs in features
        )

        n_outputs = ale_results[0].effects.shape[0]
        if is_regressor(estimator) and n_outputs > 1:
            if target is None:
                raise ValueError("target must be specified for multi-output regressors")
            if not 0 <= target < n_outputs:
                raise ValueError(f"target must be in [0, n_tasks], got {target}.")
            target_idx = target

        display = cls(
            ale_results=ale_results,
            features=features,
            feature_names=feature_names,
            target_idx=target_idx,
        )
        return display.plot(
            ax=ax, n_cols=n_cols, line_kw=line_kw, bar_kw=bar_kw, heatmap_kw=heatmap_kw
        )

    def plot(self, *, ax=None, n_cols=3, line_kw=None, bar_kw=None, heatmap_kw=None):
        """Plot accumulated local effects.

        Parameters
        ----------
        ax : Matplotlib axes or array-like of Matplotlib axes, default=None
            Axes to plot on, one per feature. If `None`, a figure is created.

        n_cols : int, default=3
            The maximum number of columns in the grid plot. Only active when
            `ax` is `None`.

        line_kw : dict, default=None
            Keyword arguments passed to ``matplotlib.pyplot.plot``.

        bar_kw : dict, default=None
            Keyword arguments passed to ``matplotlib.pyplot.bar``.

        heatmap_kw : dict, default=None
            Keyword arguments passed to ``matplotlib.pyplot.pcolormesh`` for the
            pairs of features.

        Returns
        -------
        display : :class:`~sklearn.inspection.AccumulatedLocalEffectsDisplay`
            Object that stores the computed values and the plots.
        """
        check_matplotlib_support("plot")
        import matplotlib.pyplot as plt
        from matplotlib import transforms

        line_kw = {} if line_kw is None else line_kw
        bar_kw = {} if bar_kw is None else bar_kw
        heatmap_kw = {} if heatmap_kw is None else heatmap_kw
        n_features = len(self.features)

        # Same y limits for all single-feature plots, as in PartialDependenceDisplay.
        one_way = [
            r.effects[self.target_idx] for r in self.ale_results if r.effects.ndim == 2
        ]
        if one_way:
            low, high = min(e.min() for e in one_way), max(e.max() for e in one_way)
            ylim = (low - 0.05 * (high - low), high + 0.05 * (high - low))

        if ax is None:
            n_cols = min(n_cols, n_features)
            n_rows = ceil(n_features / n_cols)
            _, axes = plt.subplots(n_rows, n_cols, squeeze=False)
            for axi in axes.ravel()[n_features:]:
                axi.set_axis_off()
        else:
            axes = np.asarray(ax, dtype=object)
            if axes.size != n_features:
                raise ValueError(
                    f"Expected ax to have {n_features} axes, got {axes.size}"
                )
            axes = axes.reshape(-1, axes.shape[-1] if axes.ndim > 1 else axes.size)

        self.axes_ = axes
        self.figure_ = axes.ravel()[0].figure
        self.lines_ = np.empty_like(axes, dtype=object)
        self.bars_ = np.empty_like(axes, dtype=object)
        self.heatmaps_ = np.empty_like(axes, dtype=object)

        for i, (fxs, result) in enumerate(zip(self.features, self.ale_results)):
            row, col = np.unravel_index(i, axes.shape)
            axi = axes[row, col]
            if len(fxs) == 2:
                self.heatmaps_[row, col] = self._plot_two_way(
                    axi, fxs, result, heatmap_kw
                )
                continue

            grid_values = result.grid_values[0]
            effects = result.effects[self.target_idx]
            if result.is_categorical[0]:
                self.bars_[row, col] = axi.bar(
                    [str(v) for v in grid_values], effects, **bar_kw
                )
                axi.tick_params(axis="x", rotation=90)
            else:
                self.lines_[row, col] = axi.plot(grid_values, effects, **line_kw)[0]
                # Ticks at the grid points, which are quantiles of the data.
                trans = transforms.blended_transform_factory(
                    axi.transData, axi.transAxes
                )
                axi.vlines(grid_values, 0, 0.05, transform=trans, color="k")
            axi.set_ylim(ylim)

            axi.set_xlabel(self.feature_names[fxs[0]])
            if col == 0:
                axi.set_ylabel("Accumulated local effect")
            else:
                axi.set_yticklabels([])

        return self

    def _plot_two_way(self, ax, features, ale_result, heatmap_kw):
        """Plot the second-order effects of a pair of features as a heatmap.

        Each cell of the grid is colored with the mean of the effects at its four
        corners. Cells without samples are masked, so they are left blank.
        """
        edges_0, edges_1 = ale_result.grid_values
        effects = ale_result.effects[self.target_idx]
        cells = (
            effects[:-1, :-1] + effects[:-1, 1:] + effects[1:, :-1] + effects[1:, 1:]
        ) / 4
        cells = np.ma.masked_where(ale_result.bin_sizes == 0, cells)

        heatmap = ax.pcolormesh(edges_0, edges_1, cells.T, **heatmap_kw)
        ax.figure.colorbar(heatmap, ax=ax, label="Accumulated local effect")
        ax.set_xlabel(self.feature_names[features[0]])
        ax.set_ylabel(self.feature_names[features[1]])
        return heatmap
