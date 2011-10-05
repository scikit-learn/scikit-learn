"""
=====================================
Blind source separation using FastICA
=====================================

:ref:`ICA` is used to estimate sources given noisy measurements.
Imagine 2 instruments playing simultaneously and 2 microphones
recording the mixed signals. ICA is used to recover the sources
ie. what is played by each instrument.

"""
print __doc__

import numpy as np
import pylab as pl
from sklearn.decomposition import FastICA

###############################################################################
# Generate sample data
np.random.seed(0)
n_samples = 2000
time = np.linspace(0, 10, n_samples)
s1 = np.sin(2 * time)  # Signal 1 : sinusoidal signal
s2 = np.sign(np.sin(3 * time))  # Signal 2 : square signal
S = np.c_[s1, s2]
S += 0.2 * np.random.normal(size=S.shape)  # Add noise

S /= S.std(axis=0)  # Standardize data
# Mix data
A = np.array([[1, 1], [0.5, 2]])  # Mixing matrix
X = np.dot(S, A.T)  # Generate observations
# Compute ICA
ica = FastICA()
S_ = ica.fit(X).transform(X)  # Get the estimated sources
A_ = ica.get_mixing_matrix()  # Get estimated mixing matrix
assert np.allclose(X, np.dot(S_, A_.T))

###############################################################################
# Plot results
pl.figure()
pl.subplot(3, 1, 1)
pl.plot(S)
pl.title('True Sources')
pl.subplot(3, 1, 2)
pl.plot(X)
pl.title('Observations (mixed signal)')
pl.subplot(3, 1, 3)
pl.plot(S_)
pl.title('ICA estimated sources')
pl.subplots_adjust(0.09, 0.04, 0.94, 0.94, 0.26, 0.36)
pl.show()
