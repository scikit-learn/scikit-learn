from __future__ import print_function

import numpy as np
import gc
from datetime import datetime
from sklearn.isotonic import IsotonicRegression
from sklearn.utils.bench import total_seconds
import sys


def bench_isotonic_regression(X, Y):
    gc.collect()

    tstart = datetime.now()
    clf = IsotonicRegression()
    clf.fit(X, Y)
    delta = datetime.now() - tstart
    return total_seconds(delta)

if __name__ == '__main__':
    min_exp, max_exp, iters = map(int, sys.argv[1:])
    timings = []
    for i in range(min_exp, max_exp):
        n = 10 ** i
        X = np.arange(n)
        Y = np.random.randint(-50, 50, size=(n,)) \
            + 50. * np.log(1 + np.arange(n))
        times = [bench_isotonic_regression(X, Y) for i in range(iters)]
        timing = (n, np.mean(times))
        timings.append(timing)
        print(n, np.mean(times))
