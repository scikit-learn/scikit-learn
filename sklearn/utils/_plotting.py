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

MULTI_PARAM_ERROR_MSG = (
    "When '{param}' is provided, it must have the same length as "
    "the number of curves to be plotted. Got: {len_param}; "
    "expected: {n_curves}."
)


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

        # We are changing from `estimator_name` to `name`, Display objects will
        # have one or the other. Try old attr name: `estimator_name` first.
        if name is None:
            name = getattr(self, "estimator_name", None)
        if name is None:
            name = getattr(self, "name", None)
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
        cls, cv_results, X, y, *, sample_weight=None, pos_label=None, fold_names=None
    ):
        check_matplotlib_support(f"{cls.__name__}.from_predictions")

        required_keys = {"estimator", "indices"}
        if not all(key in cv_results for key in required_keys):
            raise ValueError(
                "'cv_results' does not contain one of the following required keys: "
                f"{required_keys}. Set explicitly the parameters return_estimator=True "
                "and return_indices=True to the function cross_validate."
            )

        train_size, test_size = (
            len(cv_results["indices"]["train"][0]),
            len(cv_results["indices"]["test"][0]),
        )

        if _num_samples(X) != train_size + test_size:
            raise ValueError(
                "'X' does not contain the correct number of samples. "
                f"Expected {train_size + test_size}, got {_num_samples(X)}."
            )

        if type_of_target(y) != "binary":
            raise ValueError(
                f"The target y is not binary. Got {type_of_target(y)} type of"
                " target."
            )
        check_consistent_length(X, y, sample_weight)

        try:
            pos_label = _check_pos_label_consistency(pos_label, y)
        except ValueError as e:
            # Alter error message
            raise ValueError(str(e).replace("y_true", "y"))

        n_curves = len(cv_results["estimator"])
        if fold_names is None:
            fold_names = [f"Fold {idx}" for idx in range(n_curves)]
        elif len(fold_names) != n_curves:
            raise ValueError(
                MULTI_PARAM_ERROR_MSG.format(
                    param="fold_names", len_param=len(fold_names), n_curves=n_curves
                )
            )
        else:
            fold_names = fold_names

        return pos_label, fold_names

    @classmethod
    def _get_line_kwargs(
        cls,
        n_curves,
        names,
        summary_values,
        summary_value_name,
        fold_line_kwargs,
        default_line_kwargs=None,
        **kwargs,
    ):
        """Get validated line kwargs for each curve.

        Parameters
        ----------
        n_curves : int
            Number of curves.

        names : list[str]
            Names of each curve.

        summary_values : list[float]
            List of summary values for each curve (e.g., ROC AUC, average precision).

        summary_value_name : str
            Name of the summary value provided in `summary_values`.

        fold_line_kwargs : dict or list of dict
            Dictionary with keywords passed to the matplotlib's `plot` function
            to draw the individual ROC curves. If a list is provided, the
            parameters are applied to the ROC curves sequentially. If a single
            dictionary is provided, the same parameters are applied to all ROC
            curves. Ignored for single curve plots - pass as `**kwargs` for
            single curve plots.

        default_line_kwargs : dict, default=None
            Default line kwargs to be used in all curves, unless overridden by
            `fold_line_kwargs`.

        **kwargs : dict
            For a single curve plots only, keyword arguments to be passed to
            matplotlib's `plot`. Ignored for multi-curve plots - use `fold_line_kwargs`
            for multi-curve plots.
        """
        # Ensure parameters are of the correct length
        names_ = [None] * n_curves if names is None else names
        summary_values_ = (
            [None] * n_curves if summary_values is None else summary_values
        )
        # `fold_line_kwargs` ignored for single curve plots
        # `kwargs` ignored for multi-curve plots
        if n_curves == 1:
            fold_line_kwargs = [kwargs]
        else:
            fold_line_kwargs = _validate_line_kwargs(n_curves, fold_line_kwargs)

        if default_line_kwargs is None:
            default_line_kwargs = {}
        line_kwargs = []
        for fold_idx, (curve_summary_value, curve_name) in enumerate(
            zip(summary_values_, names_)
        ):
            if curve_summary_value is not None and curve_name is not None:
                default_line_kwargs["label"] = (
                    f"{curve_name} ({summary_value_name} = {curve_summary_value:0.2f})"
                )
            elif curve_summary_value is not None:
                default_line_kwargs["label"] = (
                    f"{summary_value_name} = {curve_summary_value:0.2f}"
                )
            elif curve_name is not None:
                default_line_kwargs["label"] = curve_name

            line_kwargs.append(
                _validate_style_kwargs(default_line_kwargs, fold_line_kwargs[fold_idx])
            )
        return line_kwargs


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


def _deprecate_estimator_name(old, new, version):
    """Deprecate `estimator_name` in favour of `name`."""
    version = parse_version(version)
    # Not sure if I should hard code this because this wouldn't work if we release
    # a new major version ?
    version_remove = f"{version.major}.{version.minor + 2}"
    if old != "deprecated":
        if new:
            raise ValueError(
                f"Both 'estimator_name' and 'name' provided, please only use 'name' "
                f"as 'estimator_name' is deprecated in {version} and will be removed "
                f"in {version_remove}."
            )
        warnings.warn(
            f"'estimator_name' was passed to 'name' as 'estimator_name' is deprecated "
            f"in {version} and will be removed in {version_remove}. Please use "
            f"'name' in future.",
            FutureWarning,
        )
        return old
    return new


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
        required_formatted = ", ".join(f"'{key}'" for key in required.keys())
        optional_formatted = ", ".join(f"'{key}'" for key in optional_provided.keys())
        lengths_formatted = ", ".join(
            f"{key}: {len(value)}" for key, value in all_params.items()
        )
        raise ValueError(
            f"{required_formatted}, and optional parameters {optional_formatted} "
            f"from `{class_name}` initialization (or `plot`) should all be lists of "
            f"the same length. Got: {lengths_formatted}"
        )


# Potentially useful for non binary displays `LearningCurveDisplay` and
# `ValidationCurveDisplay`, so not placed under `_BinaryClassifierCurveDisplayMixin`
def _validate_line_kwargs(n_curves, fold_line_kwargs=None, default_line_kwargs=None):
    """Ensure that `fold_names` and `fold_line_kwargs` are of correct length."""
    if fold_line_kwargs is None and default_line_kwargs is not None:
        fold_line_kwargs = default_line_kwargs
    elif fold_line_kwargs is None:
        fold_line_kwargs = [{}] * n_curves
    elif isinstance(fold_line_kwargs, Mapping):
        fold_line_kwargs = [fold_line_kwargs] * n_curves
    elif len(fold_line_kwargs) != n_curves:
        raise ValueError(
            MULTI_PARAM_ERROR_MSG.format(
                param="fold_line_kwargs",
                len_param=len(fold_line_kwargs),
                n_curves=n_curves,
            )
        )
    else:
        fold_line_kwargs = fold_line_kwargs

    return fold_line_kwargs
