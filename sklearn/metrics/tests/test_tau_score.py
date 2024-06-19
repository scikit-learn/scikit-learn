import numpy as np
from sklearn.metrics import tau_score
import pytest

def test_tau_score_perfect_prediction():
    y_true = np.array([0, 1, 0, 1])
    y_pred = np.array([0, 1, 0, 1])
    actual = tau_score(y_true, y_pred)
    expected = 1.0
    print(f"Testing perfect prediction - Expected: {expected}, Actual: {actual}")
    assert np.isclose(actual, expected), f"Expected {expected}, got {actual}"

def test_tau_score_imperfect_prediction():
    y_true = np.array([0, 1, 0, 1])
    y_pred = np.array([1, 1, 0, 0])
    actual = tau_score(y_true, y_pred)
    expected = 0.5  # Adjusted for your scenario
    print(f"Testing imperfect prediction - Expected: {expected}, Actual: {actual}")
    assert np.isclose(actual, expected), f"Expected {expected}, got {actual}"

def test_tau_score_all_wrong():
    y_true = np.array([0, 1, 2, 2])
    y_pred = np.array([1, 0, 0, 1])
    actual = tau_score(y_true, y_pred)
    expected = 0.0
    print(f"Testing all wrong predictions - Expected: {expected}, Actual: {actual}")
    assert np.isclose(actual, expected), f"Expected {expected}, got {actual}"

def test_tau_score_empty_input():
    y_true = np.array([])
    y_pred = np.array([])
    print("Testing empty input")
    with pytest.raises(ValueError):
        tau_score(y_true, y_pred)

def test_tau_score_multiclass_imperfect_prediction():
    # Define multi-class labels and imperfect predictions
    y_true = np.array([0, 1, 2, 0, 1, 2])
    y_pred = np.array([0, 2, 1, 0, 1, 1])
    
    # Calculate actual score using tau_score function
    actual = tau_score(y_true, y_pred)
    
    # Define the expected score based on previous observations or manual calculation
    expected = 0.35
    
    # Print results for debug purposes
    print(f"Testing multi-class imperfect prediction - Expected: {expected}, Actual: {actual}")
    
    # Assert to check if the actual score is close to the expected score
    assert np.isclose(actual, expected, atol=0.01), f"Expected {expected}, got {actual}"


def test_tau_score_invalid_input_type():
    y_true = ['a', 'b', 'c']
    y_pred = ['a', 'c', 'b']
    print("Testing invalid input types")
    with pytest.raises(ValueError):
        tau_score(y_true, y_pred)

if __name__ == "__main__":
    test_tau_score_perfect_prediction()
    test_tau_score_imperfect_prediction()
    test_tau_score_all_wrong()
    test_tau_score_empty_input()
    test_tau_score_invalid_input_type()
    test_tau_score_multiclass_imperfect_prediction()
