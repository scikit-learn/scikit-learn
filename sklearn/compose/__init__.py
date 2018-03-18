"""Meta-estimators for building composite models with transformers

In addition to its current contents, this module will eventually be home to
refurbished versions of Pipeline and FeatureUnion.

"""

from ._target import TransformedTargetRegressor

__all__ = [
    'TransformedTargetRegressor',
]
