# -*- coding: utf-8 -*-

import numpy as np
import NCA
import pylab as pl

# Initialisation
N = 300
aux = (np.concatenate([0.5*np.ones((N/2, 1)),
                       np.zeros((N/2, 1)), 1.1*np.ones((N/2, 1))], axis=1))
X = np.concatenate([np.random.rand(N/2, 3),
                    np.random.rand(N/2, 3) + aux])

y = np.concatenate([np.concatenate([np.ones((N/2, 1)), np.zeros((N/2, 1))]),
                    np.concatenate([np.zeros((N/2, 1)), np.ones((N/2, 1))])],
                   axis=1)
X = X.T
y = y[:, 0]
A = np.array([[1, 0, 0], [0, 1, 0]])

# Training
nca = NCA.NCA(metric=A, method='BFGS', objective='KL-divergence', options={'maxiter': 10, 'disp': True})
print nca.score(X, y)
nca.fit(X, y)
print nca.score(X, y)

# Plot
pl.subplot(2, 1, 1)
AX = np.dot(A, X)
pl.scatter(AX[0, :], AX[1, :], 30, c=y, cmap=pl.cm.Paired)
pl.subplot(2, 1, 2)
BX = np.dot(np.reshape(nca.metric, np.shape(A)), X)
pl.scatter(BX[0, :], BX[1, :], 30, c=y, cmap=pl.cm.Paired)
pl.show()
