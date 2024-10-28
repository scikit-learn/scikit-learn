"""Meta-estimators for building composite models with transformers.

In addition to its current contents, this module will eventually be home to
refurbished versions of :class:`~sklearn.pipeline.Pipeline` and
:class:`~sklearn.pipeline.FeatureUnion`.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from ._column_transformer import (
    ColumnTransformer,
    make_column_selector,
    make_column_transformer,
)
from ._target import TransformedTargetRegressor

__all__ = [
    "ColumnTransformer",
    "make_column_transformer",
    "TransformedTargetRegressor",
    "make_column_selector",
]
