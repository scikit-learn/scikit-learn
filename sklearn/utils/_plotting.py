# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause
import warnings
from collections.abc import Mapping

import numpy as np

from . import check_consistent_length
from ._optional_dependencies import check_matplotlib_support
from ._response import _get_response_values_binary
from .multiclass import type_of_target
from .validation import _check_pos_label_consistency


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

        if name is None:
            for attr in ["estimator_name", "names"]:
                name = getattr(self, attr, None)
                if name is not None:
                    break
        # One line shorter alternative, but looks funny:
        # if name is None:
        #     name = getattr(self, "estimator_name", None)
        # if name is None:
        #     name = getattr(self, "names", None)
        return ax, ax.figure, name

    @classmethod
    def _get_line_kwargs(
        cls,
        n_curves,
        names,
        summary_values,
        summary_value_name,
        fold_line_kws,
        default_line_kwargs={},
        **kwargs,
    ):
        """Get validated line kwargs for each curve."""
        # Ensure parameters are of the correct length
        names_ = [None] * n_curves if names is None else names
        summary_values_ = (
            [None] * n_curves if summary_values is None else summary_values
        )
        # `fold_line_kws` ignored for single curve plots
        # `kwargs` ignored for multi-curve plots
        if n_curves == 1:
            fold_line_kws = [kwargs]
        else:
            if fold_line_kws is None:
                # We should not set color to be the same, otherwise legend is
                # meaningless
                fold_line_kws = [
                    {"alpha": 0.5, "color": "tab:blue", "linestyle": "--"}
                ] * n_curves
            elif isinstance(fold_line_kws, Mapping):
                fold_line_kws = [fold_line_kws] * n_curves
            elif len(fold_line_kws) != n_curves:
                raise ValueError(
                    "When `fold_line_kws` is a list, it must have the same length as "
                    "the number of curves to be plotted."
                )

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
                _validate_style_kwargs(default_line_kwargs, fold_line_kws[fold_idx])
            )
        return line_kwargs

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


# TODO(1.9): remove
# Should this be a parent class method?
def _deprecate_singular(singular, plural, name):
    """Deprecate the singular version of Display parameters.

    If only `singular` parameter passed, it will be returned as a list with a warning.
    """
    if singular != "deprecated":
        warnings.warn(
            f"`{name}` was passed to `{name}s` in a list because `{name}` is "
            f"deprecated in 1.7 and will be removed in 1.9. Use "
            f"`{name}s` instead.",
            FutureWarning,
        )
        if plural:
            raise ValueError(
                f"Cannot use both `{name}` and `{name}s`. Use only `{name}s` as "
                f"`{name}` is deprecated."
            )
        return [singular]
    return plural


# Should this be a parent class/mixin method?
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
            f"from `{class_name}` initialization, should all be lists of the same "
            f"length. Got: {lengths_formatted}"
        )


def _process_fold_names_line_kwargs(n_curves, fold_names, fold_line_kwargs):
    """Ensure that `fold_names` and `fold_line_kwargs` are of correct length."""
    msg = (
        "When `{param}` is provided, it must have the same length as "
        "the number of curves to be plotted. Got {len_param} "
        "instead of {n_curves}."
    )

    if fold_names is None:
        # "<estimator> fold <idx> ?"
        fold_names_ = [f"Fold: {idx}" for idx in range(n_curves)]
    elif len(fold_names) != n_curves:
        raise ValueError(
            msg.format(param="fold_names", len_param=len(fold_names), n_curves=n_curves)
        )
    else:
        fold_names_ = fold_names

    if isinstance(fold_line_kwargs, Mapping):
        fold_line_kws_ = [fold_line_kwargs] * n_curves
    elif fold_names is not None and len(fold_line_kwargs) != n_curves:
        raise ValueError(
            msg.format(
                param="fold_line_kwargs",
                len_param=len(fold_line_kwargs),
                n_curves=n_curves,
            )
        )
    else:
        fold_line_kws_ = fold_line_kwargs

    return fold_names_, fold_line_kws_
