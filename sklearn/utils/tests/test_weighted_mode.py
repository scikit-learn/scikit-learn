import numpy as np
from sklearn.utils.extmath import weighted_mode

from scipy import stats


def test_uniform_weights():
    # with uniform weights, results should be identical to stats.mode
    x = np.random.randint(10, size=(10, 5))
    weights = np.ones(x.shape)

    for axis in (None, 0, 1):
        mode, score = stats.mode(x, axis)
        mode2, score2 = weighted_mode(x, weights, axis)

        assert np.all(mode == mode2)
        assert np.all(score == score2)


def test_random_weights():
    # set this up so that each row should have a weighted mode of 6,
    # with a score that is easily reproduced
    mode_result = 6

    x = np.random.randint(mode_result, size=(100, 10))
    w = np.random.random(x.shape)

    x[:, :5] = mode_result
    w[:, :5] += 1

    mode, score = weighted_mode(x, w, axis=1)

    assert np.all(mode == mode_result)
    assert np.all(score.ravel() == w[:, :5].sum(1))


if __name__ == '__main__':
    import nose
    nose.runmodule()
