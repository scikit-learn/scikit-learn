"""
=====================================
Bling source separation using FastICA
=====================================
"""

import numpy as np
import pylab as pl

from scikits.learn.fastica import FastICA

###############################################################################
# Generate sample data

np.random.seed(0)

n_samples = 1000
time = np.linspace(0, 10, n_samples)
s1 = np.sin(time) # Signal 1 : sinusoidal signal
s2 = np.sign(np.sin(3*time)) # Signal 2 : square signal

# Add noise

S = np.c_[s1,s2].T
S += 0.1*np.random.normal(size=S.shape)

# Standardize data
S /= S.std(axis=1)[:,np.newaxis]

# Mix data
A = [[1, 1], [0, 2]] # Mixing matrix
A = [[1, 1], [0, 2], [3, 5]] # Mixing matrix

X = np.dot(A, S) # Generate observations

# Compute ICA to reestimate sources
ica = FastICA(n_comp=2)
S_ = ica.fit(X).transform(X)

###############################################################################
# Plot results

pl.figure()

pl.subplot(3, 1, 1)
pl.plot(S.T)
pl.title('True Sources')

pl.subplot(3, 1, 2)
pl.plot(X.T)
pl.title('Observations (mixed signal)')

pl.subplot(3, 1, 3)
pl.plot(S_.T)
pl.title('ICA estimated sources')

pl.subplots_adjust(0.09, 0.04, 0.94, 0.94, 0.26, 0.36)

pl.show()