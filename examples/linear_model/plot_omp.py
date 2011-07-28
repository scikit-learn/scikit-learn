"""
===========================
Orthogonal Matching Pursuit
===========================

Using orthogonal matching pursuit for recovering a sparse signal from a noisy
measurement encoded with a dictionary
"""
print __doc__

import pylab as pl
import numpy as np
from scikits.learn.linear_model import orthogonal_mp
from scikits.learn.datasets import generate_sparse_coded_signal

n_components, n_features = 512, 100
n_atoms = 17

# generate the data
###################
y, x, D = generate_sparse_coded_signal(n_samples=1,
                                           n_components=n_components,
                                           n_features=n_features,
                                           n_nonzero_coefs=n_atoms,
                                           random_state=0)

idx, = x.nonzero()

# distort the clean signal
##########################
y_noisy = y + np.random.normal(0, 0.1, y.shape)

# plot the sparse signal
########################
pl.subplot(3, 1, 1)
pl.title("Sparse signal")
pl.stem(idx, x[idx])

# plot the noise-free reconstruction
####################################
x_r = orthogonal_mp(D.T, y, n_atoms)
idx_r, = x_r.nonzero()
pl.subplot(3, 1, 2)
pl.title("Recovered signal from noise-free measurements")
pl.stem(idx_r, x_r[idx_r])

# plot the noisy reconstruction
###############################
x_r = orthogonal_mp(D.T, y_noisy, n_atoms)
idx_r, = x_r.nonzero()
pl.subplot(3, 1, 3)
pl.title("Recovered signal from noisy measurements")
pl.stem(idx_r, x_r[idx_r])

pl.subplots_adjust(0.06, 0.04, 0.94, 0.90, 0.20, 0.38)
pl.suptitle('Sparse signal recovery with Orthogonal Matching Pursuit',
            fontsize=16)
pl.show()
