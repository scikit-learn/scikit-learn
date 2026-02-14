"""
The :mod:`sklearn._loss` module includes loss function classes suitable for
fitting classification and regression tasks.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn._loss.loss import (
    AbsoluteError,
    HalfBinomialLoss,
    HalfBinomialLossArrayAPI,
    HalfGammaLoss,
    HalfMultinomialLoss,
    HalfMultinomialLossArrayAPI,
    HalfPoissonLoss,
    HalfSquaredError,
    HalfTweedieLoss,
    HalfTweedieLossIdentity,
    HuberLoss,
    PinballLoss,
)

__all__ = [
    "AbsoluteError",
    "HalfBinomialLoss",
    "HalfBinomialLossArrayAPI",
    "HalfGammaLoss",
    "HalfMultinomialLoss",
    "HalfMultinomialLossArrayAPI",
    "HalfPoissonLoss",
    "HalfSquaredError",
    "HalfTweedieLoss",
    "HalfTweedieLossIdentity",
    "HuberLoss",
    "PinballLoss",
]
