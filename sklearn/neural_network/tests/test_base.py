import numpy as np
import pytest

from sklearn.neural_network._base import binary_log_loss, log_loss, squared_loss


def test_binary_log_loss_1_prob_finite():
    # y_proba is equal to one should result in a finite logloss
    y_true = np.array([[0, 0, 1]]).T
    y_prob = np.array([[0.9, 1.0, 1.0]]).T

    loss = binary_log_loss(y_true, y_prob, sample_weight=np.ones(shape=y_prob.shape[0]))
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
    loss = log_loss(y_true, y_prob, sample_weight=np.ones(shape=y_prob.shape[0]))
    assert np.isfinite(loss)


def test_squared_loss():
    y_true = np.array([[1, 0, 0], [0, 1, 0]])
    y_pred = np.array([[0.0, 1.0, 0.0], [0.9, 0.05, 0.05]])

    assert squared_loss(y_true, y_pred) == pytest.approx(0.3095833333333333)

    sample_weight = np.array([1, 2])
    assert squared_loss(y_true, y_pred, sample_weight) == pytest.approx(0.4525)


def test_log_loss():
    y_true = np.array([[1, 0, 0], [0, 1, 0]])
    y_pred = np.array([[0.0, 1.0, 0.0], [0.9, 0.05, 0.05]])

    assert log_loss(y_true, y_pred) == pytest.approx(19.519692831335572)

    sample_weight = np.array([1, 2])
    assert log_loss(y_true, y_pred, sample_weight) == pytest.approx(21.017558968112567)


def test_binary_log_loss():
    y_true = np.array([[0], [1]])
    y_pred = np.array([[0.1], [0.9]])

    assert (
        log_loss(y_true, y_pred)
        == pytest.approx(binary_log_loss(y_true, y_pred))
        == pytest.approx(0.10536051565782628)
    )

    sample_weight = np.array([1, 2])
    assert (
        log_loss(y_true, y_pred, sample_weight)
        == pytest.approx(binary_log_loss(y_true, y_pred, sample_weight))
        == pytest.approx(0.15804077348673942)
    )


def test_binary_log_loss_multi_label():
    y_true = np.array([[1, 0, 0], [0, 1, 0]])
    y_pred = np.array([[0.0, 1.0, 0.0], [0.9, 0.05, 0.05]])

    assert binary_log_loss(y_true, y_pred) == pytest.approx(38.71845871958495)

    sample_weight = np.array([1, 2])
    assert binary_log_loss(y_true, y_pred, sample_weight) == pytest.approx(
        41.39326405005274
    )


@pytest.mark.parametrize("loss", [squared_loss, log_loss, binary_log_loss])
def test_sample_weight_effects(loss):
    y_true = np.array([[1, 0, 0], [0, 1, 0]])
    y_pred = np.array([[0.0, 1.0, 0.0], [0.9, 0.05, 0.05]])

    sample_weight = np.array([0, 0])
    assert loss(y_true, y_pred, sample_weight) == pytest.approx(0.0)

    sample_weight = np.array([0, 1])
    assert loss(y_true, y_pred, sample_weight) == pytest.approx(
        loss(y_true[1:], y_pred[1:], sample_weight[1:])
    )

    sample_weight = np.array([1, 1])
    assert loss(y_true, y_pred, sample_weight) == pytest.approx(loss(y_true, y_pred))
