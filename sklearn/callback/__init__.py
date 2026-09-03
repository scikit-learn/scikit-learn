"""
The :mod:`sklearn.callback` module implements the framework and off the shelf
callbacks for scikit-learn estimators.
"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

from sklearn.callback._base import AutoPropagatedCallback, FitCallback
from sklearn.callback._callback_context import CallbackContext
from sklearn.callback._callback_support import CallbackSupportMixin, with_callbacks
from sklearn.callback._progressbar import ProgressBar
from sklearn.callback._scoring_monitor import ScoringMonitor, ScoringMonitorLog
from sklearn.callback.tests._utils import check_callback_estimator

__all__ = [
    "AutoPropagatedCallback",
    "CallbackContext",
    "CallbackSupportMixin",
    "FitCallback",
    "ProgressBar",
    "ScoringMonitor",
    "ScoringMonitorLog",
    "check_callback_estimator",
    "with_callbacks",
]
