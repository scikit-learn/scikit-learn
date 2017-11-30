"""
============================================================
Using PowerTransformer to transform a LogNormal distribution
============================================================

This example demonstrates the use of the Box-Cox transform
through PowerTransformer to map data from various
distributions to a Gaussian distribution.

Box-Cox is useful as a transformation in modeling problems
where homoscedasticity and normality is desired. Below
are examples of Box-Cox applied to three distributions: a
heavily left-skewed lognormal distribution, a left-skewed
chi-squared distribution, and a right-skewed Weibull
distribution.

Note that when Box-Cox is applied to the lognormal distribution,
the estimated parameter, lambda, takes on a value close to zero
since Box-Cox is defined as the natural logarithm when lambda = 0.
In all cases, the maximum likelihood estimated lambdas achieve very
close transformations to the normal distribution.
"""
print(__doc__)

# Author: Eric Chang <ericchang2017@u.northwestern.edu>
# License: BSD 3 clause

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from sklearn.preprocessing import PowerTransformer

N_SAMPLES = 3000
SEED = 42

pt = PowerTransformer(method='box-cox')
rng = np.random.RandomState(SEED)
size = (N_SAMPLES, 1)

# lognormal distribution
X_lognormal = rng.lognormal(size=size)
X_lognormal_trans = pt.fit_transform(X_lognormal)
lmbda_lognormal = round(pt.lambdas_[0], 3)

# chi-squared distribution
df = 3
X_chisq = rng.chisquare(df=df, size=size)
X_chisq_trans = pt.fit_transform(X_chisq)
lmbda_chisq = round(pt.lambdas_[0], 3)

# weibull distribution
a = 50
X_weibull = rng.weibull(a=a, size=size)
X_weibull_trans = pt.fit_transform(X_weibull)
lmbda_weibull = round(pt.lambdas_[0], 3)

# plot
params = {'font.size': 6, 'hist.bins': 150}
matplotlib.rcParams.update(params)

fig, axes = plt.subplots(nrows=2, ncols=3)
ax0, ax2, ax4, ax1, ax3, ax5 = axes.flatten()

color = 'firebrick'
ax0.hist(X_lognormal, color=color)
ax0.set_title('Lognormal')
ax1.hist(X_lognormal_trans, color=color)
title = 'Lognormal after Box-Cox, lamdba = {l}'.format(l=lmbda_lognormal)
ax1.set_title(title)

color = 'seagreen'
ax2.hist(X_chisq, color=color)
ax2.set_title('Chi-squared, df = {df}'.format(df=df))
ax3.hist(X_chisq_trans, color=color)
ax3.set_title('Chi-squared after Box-Cox, lambda = {l}'.format(l=lmbda_chisq))

color = 'royalblue'
ax4.hist(X_weibull, color=color)
ax4.set_title('Weibull, a = {a}'.format(a=a))
ax5.hist(X_weibull_trans, color=color)
ax5.set_title('Weibull after Box-Cox, lambda = {l}'.format(l=lmbda_weibull))

plt.tight_layout()
plt.show()
