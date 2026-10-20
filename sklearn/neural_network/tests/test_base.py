import numpy as np
import pytest

from sklearn._loss import HalfPoissonLoss
from sklearn.neural_network._base import binary_log_loss, log_loss, poisson_loss


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


def test_poisson_loss(global_random_seed):
    """Test Poisson loss against well tested HalfPoissonLoss."""
    n = 1000
    rng = np.random.default_rng(global_random_seed)
    y_true = rng.integers(low=0, high=10, size=n).astype(float)
    y_raw = rng.standard_normal(n)
    y_pred = np.exp(y_raw)
    sw = rng.uniform(low=0.1, high=10, size=n)

    assert 0 in y_true

    loss = poisson_loss(y_true=y_true, y_pred=y_pred, sample_weight=sw)
    pl = HalfPoissonLoss()
    loss_ref = (
        pl(y_true=y_true, raw_prediction=y_raw, sample_weight=sw)
        + pl.constant_to_optimal_zero(y_true=y_true, sample_weight=sw).mean()
        / sw.mean()
    )

    assert loss == pytest.approx(loss_ref, rel=1e-12)
