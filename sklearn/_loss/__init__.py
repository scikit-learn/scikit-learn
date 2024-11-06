"""
The :mod:`sklearn._loss` module includes loss function classes suitable for
fitting classification and regression tasks.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from .loss import (
    AbsoluteError,
    HalfBinomialLoss,
    HalfGammaLoss,
    HalfMultinomialLoss,
    HalfPoissonLoss,
    HalfSquaredError,
    HalfTweedieLoss,
    HalfTweedieLossIdentity,
    HuberLoss,
    PinballLoss,
)

__all__ = [
    "HalfSquaredError",
    "AbsoluteError",
    "PinballLoss",
    "HuberLoss",
    "HalfPoissonLoss",
    "HalfGammaLoss",
    "HalfTweedieLoss",
    "HalfTweedieLossIdentity",
    "HalfBinomialLoss",
    "HalfMultinomialLoss",
]
