"""
Benchmark for the LARS algorithm.

Work in progress
"""

from datetime import datetime
import numpy as np
from scikits.learn import glm

n, m = 100, 50000

X = np.random.randn(n, m)
y = np.random.randn(n)

if __name__ == '__main__':
    print "Computing regularization path using the LARS ..."
    start = datetime.now()
    alphas, active, path = glm.lars_path(X, y, method='lasso')
    print "This took ", datetime.now() - start


