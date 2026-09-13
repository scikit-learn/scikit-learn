"""Accumulated Local Effects (ALE) visualization class."""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from math import ceil
import numpy as np

from sklearn.inspection._accumulated_local_effects import accumulated_local_effects
from sklearn.inspection._pd_utils import _check_feature_names, _get_feature_index
from sklearn.utils._optional_dependencies import check_matplotlib_support
from sklearn.utils._plotting import _validate_style_kwargs


class AccumulatedLocalEffectsDisplay:
    """Accumulated Local Effects (ALE) Display.

    It is recommended to use
    :func:`~sklearn.inspection.AccumulatedLocalEffectsDisplay.from_estimator` to create an
    :class:`~sklearn.inspection.AccumulatedLocalEffectsDisplay`.

    Parameters
    ----------
    ale_results : Bunch
        Results of :func:`~sklearn.inspection.accumulated_local_effects`.
    features : list of int
        Indices of features for each plot.
    feature_names : list of str
        Feature names corresponding to the indices in `features`.
    target_idx : int, default=0
        In a multiclass or multioutput setting, specifies the class or target index
        for which the ALE should be plotted.

    Attributes
    ----------
    lines_ : ndarray of shape (n_plots,)
        The Line2D objects for each feature plot.
    axes_ : ndarray of shape (n_rows, n_cols)
        The matplotlib axes for the subplots.
    figure_ : matplotlib Figure
        Figure containing the plots.
    """

    def __init__(
        self,
        ale_results,
        features,
        feature_names,
        target_idx=0,
    ):
        self.ale_results = ale_results
        self.features = features
        self.feature_names = feature_names
        self.target_idx = target_idx

    def plot(
        self,
        ax=None,
        *,
        n_cols=3,
        line_kw=None,
    ):
        """Plot Accumulated Local Effects.

        Parameters
        ----------
        ax : Matplotlib axes or array-like of axes, default=None
            Axes to plot into. If None, a new figure and axes are created.
        n_cols : int, default=3
            The maximum number of columns in the grid layout.
        line_kw : dict, default=None
            Keyword arguments passed to `ax.plot`.

        Returns
        -------
        display : :class:`~sklearn.inspection.AccumulatedLocalEffectsDisplay`
            Object that stores computed display artifacts.
        """
        check_matplotlib_support(f"{self.__class__.__name__}.plot")
        import matplotlib.pyplot as plt

        n_plots = len(self.features)
        if line_kw is None:
            line_kw = {}

        if ax is None:
            n_cols = min(n_cols, n_plots)
            n_rows = int(ceil(n_plots / n_cols))
            fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 3.5 * n_rows), squeeze=False)
        elif isinstance(ax, np.ndarray):
            axes = ax.reshape(-1, n_cols) if ax.ndim == 1 else ax
            fig = axes.flat[0].figure
        else:
            axes = np.array([[ax]])
            fig = ax.figure

        self.lines_ = np.empty(n_plots, dtype=object)

        for i, feat_idx in enumerate(self.features):
            row_idx = i // n_cols
            col_idx = i % n_cols
            curr_ax = axes[row_idx, col_idx]

            x_vals = self.ale_results.grid_values[i]
            y_ale = self.ale_results.ale[i][self.target_idx]

            (line,) = curr_ax.plot(x_vals, y_ale, **line_kw)
            self.lines_[i] = line

            curr_ax.set_xlabel(self.feature_names[i])
            if col_idx == 0:
                curr_ax.set_ylabel("Accumulated Local Effect")
            curr_ax.grid(True, linestyle="--", alpha=0.6)

        # Hide any unused subplots
        for i in range(n_plots, axes.size):
            r = i // n_cols
            c = i % n_cols
            axes[r, c].set_visible(False)

        self.axes_ = axes
        self.figure_ = fig
        return self

    @classmethod
    def from_estimator(
        cls,
        estimator,
        X,
        features,
        *,
        sample_weight=None,
        feature_names=None,
        response_method="auto",
        grid_resolution=100,
        percentiles=(0.0, 1.0),
        target=None,
        ax=None,
        n_cols=3,
        line_kw=None,
    ):
        """Create and plot Accumulated Local Effects from an estimator.

        Parameters
        ----------
        estimator : BaseEstimator
            A fitted scikit-learn estimator.
        X : {array-like, dataframe} of shape (n_samples, n_features)
            The input data.
        features : array-like of {int, str} or int or str
            The feature(s) for which to compute ALE.
        sample_weight : array-like of shape (n_samples,), default=None
            Sample weights.
        feature_names : array-like of shape (n_features,), default=None
            Names of features.
        response_method : {'auto', 'predict_proba', 'decision_function', 'predict'}, default='auto'
            Prediction method.
        grid_resolution : int, default=100
            Number of intervals.
        percentiles : tuple of float, default=(0.0, 1.0)
            Extreme percentiles for interval range.
        target : int or str, default=None
            For multiclass models, the target class index or name.
        ax : Matplotlib axes, default=None
            Target axes.
        n_cols : int, default=3
            Number of columns in subplot grid.
        line_kw : dict, default=None
            Plotting keyword arguments.

        Returns
        -------
        display : :class:`~sklearn.inspection.AccumulatedLocalEffectsDisplay`
        """
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        feature_names_ = _check_feature_names(X, feature_names)

        if isinstance(features, (int, str)):
            features_list = [features]
        else:
            features_list = list(features)

        resolved_indices = [
            _get_feature_index(f, feature_names=feature_names_) for f in features_list
        ]
        resolved_names = [
            feature_names_[idx] if feature_names_ is not None else f"Feature {idx}"
            for idx in resolved_indices
        ]

        target_idx = 0
        if hasattr(estimator, "classes_") and np.size(estimator.classes_) > 2:
            if target is not None:
                if isinstance(target, str):
                    target_idx = list(estimator.classes_).index(target)
                else:
                    target_idx = int(target)

        ale_results = accumulated_local_effects(
            estimator=estimator,
            X=X,
            features=resolved_indices,
            sample_weight=sample_weight,
            feature_names=feature_names_,
            response_method=response_method,
            grid_resolution=grid_resolution,
            percentiles=percentiles,
        )

        disp = cls(
            ale_results=ale_results,
            features=resolved_indices,
            feature_names=resolved_names,
            target_idx=target_idx,
        )
        return disp.plot(ax=ax, n_cols=n_cols, line_kw=line_kw)
