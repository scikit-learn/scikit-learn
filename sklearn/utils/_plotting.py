# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
import warnings
from collections.abc import Mapping

import numpy as np

from . import check_consistent_length
from ._optional_dependencies import check_matplotlib_support
from ._response import _get_response_values_binary
from .fixes import parse_version
from .multiclass import type_of_target
from .validation import _check_pos_label_consistency, _num_samples


class _BinaryClassifierCurveDisplayMixin:
    """Mixin class to be used in Displays requiring a binary classifier.

    The aim of this class is to centralize some validations regarding the estimator and
    the target and gather the response of the estimator.
    """

    def _validate_plot_params(self, *, ax=None, name=None):
        check_matplotlib_support(f"{self.__class__.__name__}.plot")
        import matplotlib.pyplot as plt

        if ax is None:
            _, ax = plt.subplots()

        # Display classes are in process of changing from `estimator_name` to `name`.
        # Try old attr name: `estimator_name` first.
        if name is None:
            name = getattr(self, "estimator_name", getattr(self, "name", None))
        return ax, ax.figure, name

    @classmethod
    def _validate_and_get_response_values(
        cls, estimator, X, y, *, response_method="auto", pos_label=None, name=None
    ):
        check_matplotlib_support(f"{cls.__name__}.from_estimator")

        name = estimator.__class__.__name__ if name is None else name

        y_pred, pos_label = _get_response_values_binary(
            estimator,
            X,
            response_method=response_method,
            pos_label=pos_label,
        )

        return y_pred, pos_label, name

    @classmethod
    def _validate_from_predictions_params(
        cls, y_true, y_pred, *, sample_weight=None, pos_label=None, name=None
    ):
        check_matplotlib_support(f"{cls.__name__}.from_predictions")

        if type_of_target(y_true) != "binary":
            raise ValueError(
                f"The target y is not binary. Got {type_of_target(y_true)} type of"
                " target."
            )

        check_consistent_length(y_true, y_pred, sample_weight)
        pos_label = _check_pos_label_consistency(pos_label, y_true)

        name = name if name is not None else "Classifier"

        return pos_label, name

    @classmethod
    def _validate_from_cv_results_params(
        cls,
        cv_results,
        X,
        y,
        *,
        sample_weight,
        pos_label,
    ):
        check_matplotlib_support(f"{cls.__name__}.from_cv_results")

        required_keys = {"estimator", "indices"}
        if not all(key in cv_results for key in required_keys):
            raise ValueError(
                "`cv_results` does not contain one of the following required keys: "
                f"{required_keys}. Set explicitly the parameters "
                "`return_estimator=True` and `return_indices=True` to the function"
                "`cross_validate`."
            )

        train_size, test_size = (
            len(cv_results["indices"]["train"][0]),
            len(cv_results["indices"]["test"][0]),
        )

        if _num_samples(X) != train_size + test_size:
            raise ValueError(
                "`X` does not contain the correct number of samples. "
                f"Expected {train_size + test_size}, got {_num_samples(X)}."
            )

        if type_of_target(y) != "binary":
            raise ValueError(
                f"The target `y` is not binary. Got {type_of_target(y)} type of target."
            )
        check_consistent_length(X, y, sample_weight)

        try:
            pos_label = _check_pos_label_consistency(pos_label, y)
        except ValueError as e:
            # Adapt error message
            raise ValueError(str(e).replace("y_true", "y"))

        return pos_label

    @staticmethod
    def _get_legend_label(curve_legend_metric, curve_name, legend_metric_name):
        """Helper to get legend label using `name` and `legend_metric`"""
        if curve_legend_metric is not None and curve_name is not None:
            label = f"{curve_name} ({legend_metric_name} = {curve_legend_metric:0.2f})"
        elif curve_legend_metric is not None:
            label = f"{legend_metric_name} = {curve_legend_metric:0.2f}"
        elif curve_name is not None:
            label = curve_name
        else:
            label = None
        return label

    @staticmethod
    def _validate_curve_kwargs(
        n_curves,
        name,
        legend_metric,
        legend_metric_name,
        curve_kwargs,
        **kwargs,
    ):
        """Get validated line kwargs for each curve.

        Parameters
        ----------
        n_curves : int
            Number of curves.

        name : list of str or None
            Name for labeling legend entries.

        legend_metric : dict
            Dictionary with "mean" and "std" keys, or "metric" key of metric
            values for each curve. If None, "label" will not contain metric values.

        legend_metric_name : str
            Name of the summary value provided in `legend_metrics`.

        curve_kwargs : dict or list of dict or None
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual curves. If a list is provided, the
            parameters are applied to the curves sequentially. If a single
            dictionary is provided, the same parameters are applied to all
            curves.

        **kwargs : dict
            Deprecated. Keyword arguments to be passed to matplotlib's `plot`.
        """
        # TODO(1.9): Remove deprecated **kwargs
        if curve_kwargs and kwargs:
            raise ValueError(
                "Cannot provide both `curve_kwargs` and `kwargs`. `**kwargs` is "
                "deprecated in 1.7 and will be removed in 1.9. Pass all matplotlib "
                "arguments to `curve_kwargs` as a dictionary."
            )
        if kwargs:
            warnings.warn(
                "`**kwargs` is deprecated and will be removed in 1.9. Pass all "
                "matplotlib arguments to `curve_kwargs` as a dictionary instead.",
                FutureWarning,
            )
            curve_kwargs = kwargs

        if isinstance(curve_kwargs, list) and len(curve_kwargs) != n_curves:
            raise ValueError(
                f"`curve_kwargs` must be None, a dictionary or a list of length "
                f"{n_curves}. Got: {curve_kwargs}."
            )

        # Ensure valid `name` and `curve_kwargs` combination.
        if (
            isinstance(name, list)
            and len(name) != 1
            and not isinstance(curve_kwargs, list)
        ):
            raise ValueError(
                "To avoid labeling individual curves that have the same appearance, "
                f"`curve_kwargs` should be a list of {n_curves} dictionaries. "
                "Alternatively, set `name` to `None` or a single string to label "
                "a single legend entry with mean ROC AUC score of all curves."
            )

        # Ensure `name` is of the correct length
        if isinstance(name, str):
            name = [name]
        if isinstance(name, list) and len(name) == 1:
            name = name * n_curves
        name = [None] * n_curves if name is None else name

        # Ensure `curve_kwargs` is of correct length
        if isinstance(curve_kwargs, Mapping):
            curve_kwargs = [curve_kwargs] * n_curves

        default_multi_curve_kwargs = {"alpha": 0.5, "linestyle": "--", "color": "blue"}
        if curve_kwargs is None:
            if n_curves > 1:
                curve_kwargs = [default_multi_curve_kwargs] * n_curves
            else:
                curve_kwargs = [{}]

        labels = []
        if "mean" in legend_metric:
            label_aggregate = _BinaryClassifierCurveDisplayMixin._get_legend_label(
                legend_metric["mean"], name[0], legend_metric_name
            )
            # Note: "std" always `None` when "mean" is `None` - no metric value added
            # to label in this case
            if legend_metric["std"] is not None:
                # Add the "+/- std" to the end (in brackets if name provided)
                if name[0] is not None:
                    label_aggregate = (
                        label_aggregate[:-1] + f" +/- {legend_metric['std']:0.2f})"
                    )
                else:
                    label_aggregate = (
                        label_aggregate + f" +/- {legend_metric['std']:0.2f}"
                    )
            # Add `label` for first curve only, set to `None` for remaining curves
            labels.extend([label_aggregate] + [None] * (n_curves - 1))
        else:
            for curve_legend_metric, curve_name in zip(legend_metric["metric"], name):
                labels.append(
                    _BinaryClassifierCurveDisplayMixin._get_legend_label(
                        curve_legend_metric, curve_name, legend_metric_name
                    )
                )

        curve_kwargs_ = [
            _validate_style_kwargs({"label": label}, curve_kwargs[fold_idx])
            for fold_idx, label in enumerate(labels)
        ]
        return curve_kwargs_


def _validate_score_name(score_name, scoring, negate_score):
    """Validate the `score_name` parameter.

    If `score_name` is provided, we just return it as-is.
    If `score_name` is `None`, we use `Score` if `negate_score` is `False` and
    `Negative score` otherwise.
    If `score_name` is a string or a callable, we infer the name. We replace `_` by
    spaces and capitalize the first letter. We remove `neg_` and replace it by
    `"Negative"` if `negate_score` is `False` or just remove it otherwise.
    """
    if score_name is not None:
        return score_name
    elif scoring is None:
        return "Negative score" if negate_score else "Score"
    else:
        score_name = scoring.__name__ if callable(scoring) else scoring
        if negate_score:
            if score_name.startswith("neg_"):
                score_name = score_name[4:]
            else:
                score_name = f"Negative {score_name}"
        elif score_name.startswith("neg_"):
            score_name = f"Negative {score_name[4:]}"
        score_name = score_name.replace("_", " ")
        return score_name.capitalize()


def _interval_max_min_ratio(data):
    """Compute the ratio between the largest and smallest inter-point distances.

    A value larger than 5 typically indicates that the parameter range would
    better be displayed with a log scale while a linear scale would be more
    suitable otherwise.
    """
    diff = np.diff(np.sort(data))
    return diff.max() / diff.min()


def _validate_style_kwargs(default_style_kwargs, user_style_kwargs):
    """Create valid style kwargs by avoiding Matplotlib alias errors.

    Matplotlib raises an error when, for example, 'color' and 'c', or 'linestyle' and
    'ls', are specified together. To avoid this, we automatically keep only the one
    specified by the user and raise an error if the user specifies both.

    Parameters
    ----------
    default_style_kwargs : dict
        The Matplotlib style kwargs used by default in the scikit-learn display.
    user_style_kwargs : dict
        The user-defined Matplotlib style kwargs.

    Returns
    -------
    valid_style_kwargs : dict
        The validated style kwargs taking into account both default and user-defined
        Matplotlib style kwargs.
    """

    invalid_to_valid_kw = {
        "ls": "linestyle",
        "c": "color",
        "ec": "edgecolor",
        "fc": "facecolor",
        "lw": "linewidth",
        "mec": "markeredgecolor",
        "mfcalt": "markerfacecoloralt",
        "ms": "markersize",
        "mew": "markeredgewidth",
        "mfc": "markerfacecolor",
        "aa": "antialiased",
        "ds": "drawstyle",
        "font": "fontproperties",
        "family": "fontfamily",
        "name": "fontname",
        "size": "fontsize",
        "stretch": "fontstretch",
        "style": "fontstyle",
        "variant": "fontvariant",
        "weight": "fontweight",
        "ha": "horizontalalignment",
        "va": "verticalalignment",
        "ma": "multialignment",
    }
    for invalid_key, valid_key in invalid_to_valid_kw.items():
        if invalid_key in user_style_kwargs and valid_key in user_style_kwargs:
            raise TypeError(
                f"Got both {invalid_key} and {valid_key}, which are aliases of one "
                "another"
            )
    valid_style_kwargs = default_style_kwargs.copy()

    for key in user_style_kwargs.keys():
        if key in invalid_to_valid_kw:
            valid_style_kwargs[invalid_to_valid_kw[key]] = user_style_kwargs[key]
        else:
            valid_style_kwargs[key] = user_style_kwargs[key]

    return valid_style_kwargs


def _despine(ax):
    """Remove the top and right spines of the plot.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes of the plot to despine.
    """
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    for s in ["bottom", "left"]:
        ax.spines[s].set_bounds(0, 1)


def _deprecate_estimator_name(estimator_name, name, version):
    """Deprecate `estimator_name` in favour of `name`."""
    version = parse_version(version)
    version_remove = f"{version.major}.{version.minor + 2}"
    if estimator_name != "deprecated":
        if name:
            raise ValueError(
                "Cannot provide both `estimator_name` and `name`. `estimator_name` "
                f"is deprecated in {version} and will be removed in {version_remove}. "
                "Use `name` only."
            )
        warnings.warn(
            f"`estimator_name` is deprecated in {version} and will be removed in "
            f"{version_remove}. Use `name` instead.",
            FutureWarning,
        )
        return estimator_name
    return name


def _convert_to_list_leaving_none(param):
    """Convert parameters to a list, leaving `None` as is."""
    if param is None:
        return None
    if isinstance(param, list):
        return param
    return [param]


def _check_param_lengths(required, optional, class_name):
    """Check required and optional parameters are of the same length."""
    optional_provided = {}
    for name, param in optional.items():
        if isinstance(param, list):
            optional_provided[name] = param

    all_params = {**required, **optional_provided}
    if len({len(param) for param in all_params.values()}) > 1:
        param_keys = [key for key in all_params.keys()]
        # Note: below code requires `len(param_keys) >= 2`, which is the case for all
        # display classes
        params_formatted = " and ".join([", ".join(param_keys[:-1]), param_keys[-1]])
        or_plot = ""
        if "'name' (or self.name)" in param_keys:
            or_plot = " (or `plot`)"
        lengths_formatted = ", ".join(
            f"{key}: {len(value)}" for key, value in all_params.items()
        )
        raise ValueError(
            f"{params_formatted} from `{class_name}` initialization{or_plot}, "
            f"should all be lists of the same length. Got: {lengths_formatted}"
        )
