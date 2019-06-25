"""Meta-estimators for building composite models with transformers

In addition to its current contents, this module will eventually be home to
refurbished versions of Pipeline and FeatureUnion.

"""

from ._column_transformer import ColumnTransformer, make_column_transformer
from ._target import TransformedTargetRegressor
from ._resampled import ResampledTrainer


__all__ = [
    'ColumnTransformer',
    'make_column_transformer',
    'TransformedTargetRegressor',
    'ResampledTrainer',
]
