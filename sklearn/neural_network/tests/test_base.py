import numpy as np
import pytest

from sklearn.neural_network._base import binary_log_loss, log_loss, logcosh_loss


def test_binary_log_loss_1_prob_finite():
    # y_proba is equal to one should result in a finite logloss
    y_true = np.array([[0, 0, 1]]).T
    y_prob = np.array([[0.9, 1.0, 1.0]]).T

    loss = binary_log_loss(y_true, y_prob)
    assert np.isfinite(loss)


@pytest.mark.parametrize(
    "y_true, y_prob",
    [
        (
            np.array([[1, 0, 0], [0, 1, 0]]),
            np.array([[0.0, 1.0, 0.0], [0.9, 0.05, 0.05]]),
        ),
        (np.array([[0, 0, 1]]).T, np.array([[0.9, 1.0, 1.0]]).T),
    ],
)
def test_log_loss_1_prob_finite(y_true, y_prob):
    # y_proba is equal to 1 should result in a finite logloss
    loss = log_loss(y_true, y_prob)
    assert np.isfinite(loss)


def test_logcosh_loss():
    """Test logcosh_loss with various error values."""

    # Test mixed positive and negative error
    y_true = np.array([1, -1, 2, -2, 0])
    y_pred = np.array([0, 0, 0, 0, 0])
    expected_loss = 0.7035134  # np.log(np.cosh(np.array([1, -1, 2, -2, 0]))).mean()
    loss = logcosh_loss(y_true, y_pred)
    assert np.allclose(loss, expected_loss, rtol=1e-04)
