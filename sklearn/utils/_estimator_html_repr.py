# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

import warnings

from ._repr_html.base import _HTMLDocumentationLinkMixin
from ._repr_html.estimator import (
    _get_visual_block,
    _IDCounter,
    _VisualBlock,
    _write_estimator_html,
    _write_label_html,
    estimator_html_repr,
)

__all__ = [
    "_HTMLDocumentationLinkMixin",
    "_IDCounter",
    "_VisualBlock",
    "_get_visual_block",
    "_write_estimator_html",
    "_write_label_html",
    "estimator_html_repr",
]

# TODO(1.8): Remove the entire module
warnings.warn(
    "Importing from sklearn.utils._estimator_html_repr is deprecated. The tools have "
    "been moved to sklearn.utils._repr_html. Be aware that this module is private and "
    "may be subject to change in the future. The module _estimator_html_repr will be "
    "removed in 1.8.0.",
    FutureWarning,
    stacklevel=2,
)
