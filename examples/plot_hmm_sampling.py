"""
==================================
Demonstration of sampling from HMM
==================================

This script shows how to sample points from HMM.
"""

import numpy as np
from sklearn import hmm
import matplotlib.pyplot as plt

##############################################################
# prepareing parameters
startprob = np.array([0.6, 0.3, 0.1])
transmat = np.array([[0.7, 0.2, 0.1], [0.3, 0.5, 0.2],
                [0.2, 0.2, 0.6]])
means = np.array([[0.0, 0.0], [5.0, -1.0], [5.0, 10.0]])
covars = np.tile(np.identity(2), (3, 1, 1))

# build an HMM instance and set parameters
model = hmm.GaussianHMM(3, "full", startprob, transmat)
model.means_ = means
model.covars_ = covars
###############################################################

# generate samples
X, Z = model.sample(500)

#plot the sampled data
plt.plot(X[:, 0], X[:, 1], "-o", label="observable", ms=10,
        mfc="orange", alpha=0.7)
plt.legend()
plt.show()
