import warnings

warnings.warn("sklearn.feature_selection.selector_mixin.SelectorMixin "
              "has been renamed "
              "sklearn.feature_selection.from_model._LearntSelectorMixin, "
              "and this alias will be removed in version 0.16",
              DeprecationWarning)

from .from_model import _LearntSelectorMixin as SelectorMixin

__all__ = ['SelectorMixin']
