import pytest
import numpy as np
from numpy.testing import assert_allclose


def test_fft_function():
    # Many NumPy symbols are imported into the scipy namespace, including
    # numpy.fft.fft as scipy.fft, conflicting with this module (gh-10253)
    np.random.seed(1234)

    # Callable before scipy.fft is imported
    import scipy
    x = np.random.randn(10) + 1j * np.random.randn(10)
    with pytest.deprecated_call(match=r'1\.5\.0'):
        X = scipy.fft(x)
    with pytest.deprecated_call(match=r'2\.0\.0'):
        y = scipy.ifft(X)
    assert_allclose(y, x)

    # Callable after scipy.fft is imported
    import scipy.fft
    assert_allclose(X, scipy.fft.fft(x))
    with pytest.deprecated_call(match=r'1\.5\.0'):
        X = scipy.fft(x)
    assert_allclose(X, scipy.fft.fft(x))
    with pytest.deprecated_call(match=r'2\.0\.0'):
        y = scipy.ifft(X)
    assert_allclose(y, x)

    # Callable when imported using from
    from scipy import fft
    with pytest.deprecated_call(match=r'1\.5\.0'):
        X = fft(x)
    with pytest.deprecated_call(match=r'2\.0\.0'):
        y = scipy.ifft(X)
    assert_allclose(y, x)
