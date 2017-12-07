import pytest
import time

from sklearn.utils.random import loguniform
from scipy.stats import expon

def timer(func, count):
    start = time.time()

    for _ in range(count):
        _ = func()

    end = time.time()
    avg = (end - start) / count

    return avg

def test_time():
    ex = expon()
    lu = loguniform()

    assert timer(ex.rvs, 1000) < timer(lu.rvs, 1000)
