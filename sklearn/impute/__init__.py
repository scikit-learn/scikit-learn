"""Transformers for missing value imputation"""
import typing
import sys

from ._base import MissingIndicator, SimpleImputer
from ._knn import KNNImputer

if typing.TYPE_CHECKING:
    # Avoid errors in type checkers (e.g. mypy) for experimental estimators.
    # TODO: remove this check once the estimator is no longer experimental.
    from ._iterative import IterativeImputer  # noqa

__all__ = ["MissingIndicator", "SimpleImputer", "KNNImputer"]


# TODO: remove this check once the estimator is no longer experimental.
def __getattr__(name):
    if name == "IterativeImputer":
        if "sklearn.experimental.enable_iterative_imputer" not in sys.modules:
            raise ImportError(
                f"{name} is experimental and the API might change without any "
                "deprecation cycle. To use it, you need to explicitly import "
                "enable_iterative_imputer:\n"
                "from sklearn.experimental import enable_iterative_imputer"
            )
        from . import _iterative

        return getattr(_iterative, name)
    raise AttributeError(f"module {__name__} has no attribute {name}")
