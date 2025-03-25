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

AGGREGATE_ERROR_MESSAGE = (
    "'fold_line_kwargs' must be a single dictionary to be applied to all curves "
    "when {param_value}, as only one legend entry will be added."
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
        name,
        fold_line_kwargs,
        show_aggregate_score,
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
        if show_aggregate_score:
            if isinstance(name, list) and len(name) != 1:
                raise ValueError(
                    "'name' must be a string or list of length one when "
                    "'show_aggregate_score' is True, as only one legend entry "
                    "will be added."
                )
            # Should we allow a list of length 1 (in addition to single dict) ??
            if isinstance(fold_line_kwargs, list):
                raise ValueError(
                    AGGREGATE_ERROR_MESSAGE.format(
                        param_value="'show_aggregate_score' is True"
                    )
                )
        else:
            # Individual ROC AUC scores shown
            if name is not None and (isinstance(name, list) and len(name) != n_curves):
                raise ValueError(
                    f"'name' must be None or list of length {n_curves} when "
                    f"'show_aggregate_score' is False."
                )
            if name is None:
                name = [f"Fold {idx}" for idx in range(n_curves)]
        return pos_label, name

    @classmethod
    def _get_legend_label(cls, curve_summary_value, curve_name, summary_value_name):
        """Helper to get legend label using `name_` and `summary_value_`"""
        if curve_summary_value is not None and curve_name is not None:
            label = f"{curve_name} ({summary_value_name} = {curve_summary_value:0.2f})"
        elif curve_summary_value is not None:
            label = f"{summary_value_name} = {curve_summary_value:0.2f}"
        elif curve_name is not None:
            label = curve_name
        return label

    @classmethod
    def _get_line_kwargs(
        cls,
        n_curves,
        name,
        summary_value,
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

        name : list[str] or None
            Name for labeling legend entries.

        summary_value : list[float] or tuple(float, float) or None
            Either list of `n_curves` summary values for each curve (e.g., ROC AUC,
            average precision) or a single float summary value for all curves.

        summary_value_name : str or None
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
        name_ = [None] * n_curves if name is None else name
        summary_value_ = [None] * n_curves if summary_value is None else summary_value
        # `fold_line_kwargs` ignored for single curve plots
        # `kwargs` ignored for multi-curve plots
        if n_curves == 1:
            fold_line_kwargs = [kwargs]
        else:
            # Should we add an extra check ensuring `fold_line_kwargs` is the
            # same when `summary_value` is (mean, std) ? or allow this flexibility
            # when using `plot` - note this is checked and prevented in
            # `from_cv_results`
            fold_line_kwargs = cls._validate_line_kwargs(n_curves, fold_line_kwargs)

        if default_line_kwargs is None:
            default_line_kwargs = {}

        labels = []
        if isinstance(summary_value_, tuple):
            label_aggregate = cls._get_legend_label(
                summary_value_[0], name_[0], summary_value_name
            )
            label_aggregate = label_aggregate + f" +/- {summary_value_[1]:0.2f}"
            labels.extend([label_aggregate] + [None] * (n_curves - 1))
        else:
            for curve_summary_value, curve_name in zip(summary_value_, name_):
                labels.append(
                    cls._get_legend_label(
                        curve_summary_value, curve_name, summary_value_name
                    )
                )

        line_kwargs = []
        for fold_idx, label in enumerate(labels):
            default_line_kwargs["label"] = label
            line_kwargs.append(
                _validate_style_kwargs(default_line_kwargs, fold_line_kwargs[fold_idx])
            )
        return line_kwargs

    # TODO: if useful for non binary displays (e.g.,`LearningCurveDisplay`,
    # `ValidationCurveDisplay`) amend to function
    @classmethod
    def _validate_line_kwargs(
        cls,
        n_curves,
        fold_line_kwargs=None,
        default_line_kwargs=None,
        aggregate_error=False,
    ):
        """Ensure `fold_line_kwargs` length and incorporate default kwargs.

        * If `fold_line_kwargs` is None:
          * If `default_line_kwargs` is None, list of `n_curves` empty dictionaries
            is returned.
          * If `default_line_kwargs` is not None, list of `n_curves` dictionaries
            of `default_line_kwargs` returned.
        * If `fold_line_kwargs` is a single dictionary, it is incorporated with
          `default_line_kwargs` using `_validate_style_kwargs`, and the resulting
          dictionary is repeated `n_curves` times and returned.
        * If `fold_line_kwargs` is a list of length `n_curves`, each dict is
          incorporated with `default_line_kwargs` using `_validate_style_kwargs` and
          returned as list of `n_curves` dictionaries.

        If `fold_line_kwargs` is a list not of length `n_curves`, an error is raised.

        Parameters
        ----------
        n_curves : int
            Number of curves.

        fold_line_kwargs : dict or list of dict, default=None
            Keywords arguments to be passed to matplotlib's `plot` function
            to draw ROC curves.

        default_line_kwargs : dict, default=None
            Default line kwargs to be used in all curves, unless overridden by
            `fold_line_kwargs`.

        Returns
        -------
        fold_line_kwargs : list of dict
            List of `n_curves` dictionaries.
        """
        if fold_line_kwargs is None and default_line_kwargs is None:
            fold_line_kwargs = [{}] * n_curves
        elif fold_line_kwargs is None and default_line_kwargs is not None:
            fold_line_kwargs = default_line_kwargs

        if isinstance(fold_line_kwargs, Mapping):
            fold_line_kwargs = [fold_line_kwargs] * n_curves
        elif len(fold_line_kwargs) != n_curves:
            raise ValueError(
                f"'fold_line_kwargs' must be a list of length {n_curves} or a "
                f"single dictionary. Got list of length: {len(fold_line_kwargs)}."
            )
        else:
            fold_line_kwargs = fold_line_kwargs

        if default_line_kwargs is not None:
            fold_line_kwargs = [
                _validate_style_kwargs(default_line_kwargs, single_kwargs)
                for single_kwargs in fold_line_kwargs
            ]

        return fold_line_kwargs


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
    version_remove = f"{version.major}.{version.minor + 2}"
    if old != "deprecated":
        if new:
            raise ValueError(
                f"Both 'estimator_name' and 'name' provided, please only use 'name' "
                f"as 'estimator_name' is deprecated in {version} and will be removed "
                f"in {version_remove}."
            )
        warnings.warn(
            f"'estimator_name' is deprecated in {version} and will be removed in "
            f"{version_remove}. The value of 'estimator_name' was passed to 'name'"
            "but please use 'name' in future.",
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
