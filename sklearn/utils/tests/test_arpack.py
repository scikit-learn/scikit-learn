import numpy as np

from sklearn.utils.arpack import _init_arpack_v0


def test_init_arpack_v0():
    v0s = [_init_arpack_v0(1000, 0)]
    for i in range(1, 100):
        v0s.append(_init_arpack_v0(1000, i))
        assert not any(np.equal(v0s[i], v0s[i-1]))

    v0 = np.concatenate(v0s)
    assert np.allclose(np.mean(v0), 0, atol=1e-2)
    assert np.allclose(np.std(v0), 1/np.sqrt(3), atol=1e-3)
