"""
The :mod:`sklearn._loss` module includes loss function classes suitable for
fitting classification and regression tasks.
"""

from .loss import (
    HalfSquaredError,
    AbsoluteError,
    PinballLoss,
    HalfPoissonLoss,
    HalfGammaLoss,
    HalfTweedieLoss,
    HalfTweedieLossIdentity,
    HalfBinomialLoss,
    HalfMultinomialLoss,
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
