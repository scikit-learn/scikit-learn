#!/usr/bin/python
# -*- coding: utf-8 -*-

r"""
============================================================================
Matern kernel: influence of coef0 on kernel covariance
============================================================================

The example shows how the kernel covariance decreases with increasing
dissimilarity of the two inputs for different values of coef0 (the parameter
"nu" of the Matern kernel)

See Rasmussen and Williams 2006, pp84 for details regarding the different
variants of the Matern kernel.

"""
print(__doc__)

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause


import numpy as np

import matplotlib.pyplot as plt

from sklearn.metrics.pairwise import matern_kernel

d = np.linspace(-4, 4, 500)[:, None]

for coef0 in [0.5, 1.5, 2.5, np.inf]:
    K = matern_kernel(d, [[0.0]], gamma=1, coef0=coef0)
    plt.plot(d[:, 0], K[:, 0], label=coef0)

plt.xlabel("distance")
plt.ylabel("covariance")
plt.yscale("log")
plt.ylim(1e-3, 1e0)
plt.legend(title="coef0")
plt.show()
