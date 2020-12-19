import pytest
import numpy as np

from sklearn.neural_network._base import binary_log_loss, squared_loss
from sklearn.neural_network._base import log_loss


def test_binary_log_loss_1_prob_finite():
    # y_proba is equal to one should result in a finite logloss
    y_true = np.array([[0, 0, 1]]).T
    y_prob = np.array([[0.9, 1.0, 1.0]]).T

    loss = binary_log_loss(y_true, y_prob,
                           sample_weight=np.ones(y_true.shape[0]))
    assert np.isfinite(loss)


@pytest.mark.parametrize("y_true, y_prob", [
    (np.array([[1, 0, 0], [0, 1, 0]]),
     np.array([[0., 1., 0.], [0.9, 0.05, 0.05]])),
    (np.array([[0, 0, 1]]).T,
     np.array([[0.9, 1.0, 1.0]]).T),
])
def test_log_loss_1_prob_finite(y_true, y_prob):
    # y_proba is equal to 1 should result in a finite logloss
    loss = log_loss(y_true, y_prob, sample_weight=np.ones(y_true.shape[0]))
    assert np.isfinite(loss)


def test_loss_functions():
    # for testing the correctness of all loss functions

    # squared loss
    y_true = np.array([[1, 2, 3]]).T
    y_pred = np.array([[1.2, 1.8, 3.6]]).T

    square_loss_value = squared_loss(y_true, y_pred,
                                     sample_weight=np.ones(y_true.shape[0]))
    exp_square_loss_value = 0.07333333333333335
    assert square_loss_value == exp_square_loss_value

    # binary log-loss
    y_true = np.array([[0, 1, 1]]).T
    y_prob = np.array([[0.4, 0.6, 0.95]]).T

    binary_loss_value = binary_log_loss(y_true, y_prob,
                                        sample_weight=np.ones(y_true.shape[0]))
    exp_binary_loss_value = 0.357648180639844
    assert binary_loss_value == exp_binary_loss_value

    # log loss
    y_true = np.array([[0, 1, 2]])
    y_prob = np.array([[0.8, 0.15, 0.05], [0.2, 0.65, 0.15], [0.2, 0.3, 0.5]])
    log_loss_value = log_loss(y_true, y_prob,
                              sample_weight=np.ones(y_true.shape[0]))
    exp_log_loss_value = 4.901291527767969
    assert log_loss_value == exp_log_loss_value
