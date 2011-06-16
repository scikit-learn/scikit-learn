"""Using Orthogonal Matching Pursuit for recovering a sparse signal from
a noisy measurement encoded with a dictionary"""

import pylab as pl
import numpy as np
from scikits.learn.linear_model import orthogonal_mp

n_samples, n_features = 512, 100
n_atoms = 17

# generate the dictionary
#########################
X = np.random.randn(n_samples, n_features)
X /= np.sqrt(np.sum((X ** 2), axis=0))

# generate the sparse signal
############################
gamma = np.zeros(n_features)
idx = np.arange(n_features)
np.random.shuffle(idx)
idx = idx[:n_atoms]
gamma[idx] = np.random.normal(0, 5, n_atoms)

# encode the signal
###################
y = np.dot(X, gamma)
y_noisy = y + np.random.normal(0, 0.1, y.shape)

# plot the sparse signal
########################
pl.subplot(3, 1, 1)
pl.title("Sparse signal")
pl.stem(idx, gamma[idx])

# plot the noise-free reconstruction
####################################
x_r = orthogonal_mp(X, y, n_atoms)
idx_r = np.where(x_r != 0)[0]
print idx_r
print x_r[idx_r]
pl.subplot(3, 1, 2)
pl.title("Recovered signal from noise-free measurements")
pl.stem(idx_r, x_r[idx_r])

# plot the noisy reconstruction
###############################
x_r = orthogonal_mp(X, y_noisy, n_atoms)
idx_r = np.where(x_r != 0)[0]
pl.subplot(3, 1, 3)
pl.title("Recovered signal from noisy measurements")
pl.stem(idx_r, x_r[idx_r])

pl.subplots_adjust(0.06, 0.04, 0.94, 0.90, 0.20, 0.38)
pl.suptitle('Sparse signal recovery with Orthogonal Matching Pursuit',
            fontsize=16)
pl.show()