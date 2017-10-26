"""
==================================================================
Finding the driving force of a nonlinear dynamical system with SFA
==================================================================

This example reproduces some of the
results reported in Laurenz Wiskott, `Estimating Driving Forces of
Nonstationary Time Series with Slow Feature Analysis`
(`arXiv.org e-Print archive <http://arxiv.org/abs/cond-mat/0312317>`_).

We record a time series generated from a non-linear dynamical system,
driven by a non-stationary but slowly varying driving force.

Slow Feature Analysis (SFA) finds a set of uncorrelated directions in input
space, along which the signal varies as slowly as possible, which makes it
possible to reconstruct the driving force, starting from the chaotic-looking
time series.

SFA is a linear algorithm, but it can be easily extended to the nonlinear case
by expanding the input with a set of nonlinear basis.
"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.decomposition import SFA


print(__doc__)

# Define the underlying driving force of the chaotic system.
# Time axis: 1 second, with a sample rate of 10 KHz
t = np.linspace(0, 1, 10000, endpoint=False)
driving_force = (
    np.sin(2 * np.pi * 5 * t) +
    np.sin(2 * np.pi * 11 * t) +
    np.sin(2 * np.pi * 13 * t)
)

# Use the driving force as the parameter to our non-linear dynamical system.
def logistic_map(x, r):
    return r * x * (1 - x)

series = np.zeros((10000, 1), dtype=np.float64)
series[0] = 0.6  # initial ocndition
for i in range(1, 10000):
    series[i] = logistic_map(series[i - 1], 3.6 + 0.13 * driving_force[i])

# Plot the input time series.
plt.figure()
plt.plot(series, 'k.', markersize=1)
plt.title('Input time series')

# The time series is measured along a sliding window, and expanded in the
# space of polynomials of degree 3.
window_size = 9
X_time_window = np.lib.stride_tricks.as_strided(
    series, shape=(series.shape[0] - window_size, window_size), strides=(8, 8))

pipeline = Pipeline([
    ('expansion', PolynomialFeatures(degree=3, include_bias=False)),
    ('sfa', SFA(n_components=1)),
])

# Find the slowest-varying signal that can be extracted from the input.
reconstructed = pipeline.fit_transform(X_time_window.copy())

normalized_force = (driving_force - driving_force.mean()) / driving_force.std()
corr = np.corrcoef(reconstructed,
                   normalized_force[:-9, np.newaxis], rowvar=False)
print('Correleation coefficient of driving force with slowest signal: '
      '{:.3f}'.format(corr[0, 1]))

plt.figure()
plt.plot(reconstructed, 'r--')
plt.plot(normalized_force, 'k-')
plt.legend(['Driving force', 'Slowest signal'])
plt.title('Comparison of the driving force and its reconstruction')

plt.show()
