"""
The :mod:`sklearn._loss` module includes loss function classes suitable for
fitting classification and regression tasks.
"""

from ._loss import (
    HalfSquaredError,
    AbsoluteError,
    PinballLoss,
    HalfPoissonLoss,
    HalfGammaLoss,
    HalfTweedieLoss,
    BinaryCrossEntropy,
    CategoricalCrossEntropy,
)


__all__ = [
    "HalfSquaredError",
    "AbsoluteError",
    "PinballLoss",
    "HalfPoissonLoss",
    "HalfGammaLoss",
    "HalfTweedieLoss",
    "BinaryCrossEntropy",
    "CategoricalCrossEntropy",
]
