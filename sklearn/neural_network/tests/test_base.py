import numpy as np
import pytest

from sklearn.neural_network._base import binary_log_loss, log_loss, squared_loss


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


def test_squared_loss():
    y_true = np.array([[0, 1], [0, 0]])
    y_pred = np.array([[1.0, 1.0], [1.0, 0.0]])

    assert squared_loss(y_true, y_pred) == pytest.approx(0.25)

    sample_weight = np.array([1, 2])
    assert squared_loss(y_true, y_pred, sample_weight) == pytest.approx(0.375)


def test_log_loss():
    y_true = np.array([[1, 0, 0], [0, 1, 0]])
    y_prob = np.exp([[-0.2, -2, -2], [-1, -0.1, -1]])

    assert log_loss(y_true, y_prob) == pytest.approx(0.15)

    sample_weight = np.array([1, 2])
    assert log_loss(y_true, y_prob, sample_weight) == pytest.approx(0.2)


def test_binary_log_loss():
    y_true = np.array([[0], [1]])
    y_prob = np.exp([[-2], [-0.2]])

    assert (
        log_loss(y_true, y_prob)
        == pytest.approx(binary_log_loss(y_true, y_prob))
        == pytest.approx(0.17270672893442957)
    )

    sample_weight = np.array([1, 2])
    assert (
        log_loss(y_true, y_prob, sample_weight)
        == pytest.approx(binary_log_loss(y_true, y_prob, sample_weight))
        == pytest.approx(0.2727067289344296)
    )


def test_binary_log_loss_multi_label():
    y_true = np.array([[1, 0, 0], [0, 1, 0]])
    y_prob = np.exp([[-0.2, -2, -2], [-1, -0.1, -1]])

    assert binary_log_loss(y_true, y_prob) == pytest.approx(0.7540886032559411)

    sample_weight = np.array([1, 2])
    assert binary_log_loss(y_true, y_prob, sample_weight) == pytest.approx(
        1.262763748643023
    )


@pytest.mark.parametrize(
    "y_true, y_prob",
    [
        (
            np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            np.array([[0.0, 1.0, 0.0], [0.9, 0.05, 0.05], [0.1, 0.1, 0.8]]),
        ),
        (np.array([[0, 0, 1]]).T, np.array([[0.9, 1.0, 1.0]]).T),
    ],
)
@pytest.mark.parametrize("loss", [squared_loss, log_loss, binary_log_loss])
def test_sample_weight_ones_effect(y_true, y_prob, loss):
    sample_weight = np.array([1, 1, 1])
    assert loss(y_true, y_prob, sample_weight) == pytest.approx(loss(y_true, y_prob))
