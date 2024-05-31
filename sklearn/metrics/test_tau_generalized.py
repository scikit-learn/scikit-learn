import numpy as np
import pytest
from tau_generalized import TauGeneralized
from sklearn.metrics import confusion_matrix

def test_perfect_prediction():
    y_true = [0, 1, 2, 0, 1, 2]
    y_pred = [0, 1, 2, 0, 1, 2]
    tau = TauGeneralized(y_true, y_pred)
    assert tau.get_tau() == 1, "Tau should be 1 for perfect prediction"

def test_random_guess():
    y_true = [0, 1, 2, 0, 1, 2]
    y_pred = [2, 1, 0, 1, 0, 2]
    tau = TauGeneralized(y_true, y_pred)
    # Expecting a Tau score lower than 1, adjust according to your model's expected performance
    assert tau.get_tau() < 1, "Tau should be less than 1 for random guesses"

def test_worst_case():
    y_true = [0, 0, 0, 1, 1, 1]
    y_pred = [1, 1, 1, 0, 0, 0]
    tau = TauGeneralized(y_true, y_pred)
    assert tau.get_tau() == 0, "Tau should be 0 for completely incorrect predictions"

def test_with_normalization():
    y_true = [0, 1, 0, 1]
    y_pred = [0, 1, 1, 0]
    tau = TauGeneralized(y_true, y_pred, normalize=True)
    # Check some expected value that makes sense after normalization
    assert 0 < tau.get_tau() < 1, "Normalized Tau should be between 0 and 1"

# You can add more tests as necessary
