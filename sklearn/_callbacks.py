# License: BSD 3 clause
from typing import List, Callable, Optional
from abc import ABC, abstractmethod

import numpy as np

CALLBACK_PARAM_TYPES = {
    "n_iter": int,
    "max_iter": int,
    "loss": (float, dict),
    "score": (float, dict),
    "validation_loss": (float, dict),
    "validation_score": (float, dict),
    "coef": np.ndarray,
    "intercept": (np.ndarray, float),
}


def _check_callback_params(**kwargs):
    invalid_params = []
    invalid_types = []
    for key, val in kwargs.items():
        if key not in CALLBACK_PARAM_TYPES:
            invalid_params.append(key)
        else:
            val_types = CALLBACK_PARAM_TYPES[key]
            if not isinstance(val, val_types):
                invalid_types.append(f"{key}={val} is not of type {val_types}")
    msg = ""
    if invalid_params:
        msg += ("Invalid callback parameters: {}, must be one of {}. ").format(
            ", ".join(invalid_params), ", ".join(CALLBACK_PARAM_TYPES.keys())
        )
    if invalid_types:
        msg += "Invalid callback parameters: " + ", ".join(invalid_types)
    if msg:
        raise ValueError(msg)


def _eval_callbacks(
    callbacks: Optional[List[Callable]], method="on_iter_end", **kwargs
) -> None:
    if callbacks is None:
        return

    for callback in callbacks:
        getattr(callback, method)(**kwargs)


class BaseCallback(ABC):
    @abstractmethod
    def on_fit_begin(self, estimator, X, y) -> None:
        pass

    @abstractmethod
    def on_iter_end(self, **kwargs) -> None:
        pass
