"""
The :mod:`sklearn._loss` module includes loss function classes suitable for
fitting classification and regression tasks.
"""

from .loss import (
    AbsoluteError,
    HalfBinomialLoss,
    HalfGammaLoss,
    HalfMultinomialLoss,
    HalfPoissonLoss,
    HalfSquaredError,
    HalfTweedieLoss,
    HalfTweedieLossIdentity,
    PinballLoss,
)

__all__ = [
    "HalfSquaredError",
    "AbsoluteError",
    "PinballLoss",
    "HalfPoissonLoss",
    "HalfGammaLoss",
    "HalfTweedieLoss",
    "HalfTweedieLossIdentity",
    "HalfBinomialLoss",
    "HalfMultinomialLoss",
]
