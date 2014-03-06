import numpy as np
import gc
from datetime import datetime
from sklearn.isotonic import IsotonicRegression
from sklearn.utils.bench import total_seconds
import sys
import matplotlib.pyplot as plt


def bench_isotonic_regression(X, Y):
    gc.collect()

    # start time
    tstart = datetime.now()
    clf = IsotonicRegression()
    clf.fit(X, Y)
    delta = (datetime.now() - tstart)
    return total_seconds(delta)

if __name__ == '__main__':
    n, iters = int(sys.argv[1]), int(sys.argv[2])
    X = np.arange(n)
    rng = np.random.RandomState(42)
    Y = rng.randint(-50, 50, size=(n,)) + 50. * np.log(1 + np.arange(n))
    times = [bench_isotonic_regression(X, Y) for i in range(iters)]
    print times
    print np.mean(times)
    plt.hist(times)
    plt.show()
