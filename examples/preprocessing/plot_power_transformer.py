"""
==========================================================
Using PowerTransformer to apply the Box-Cox transformation
==========================================================

This example demonstrates the use of the Box-Cox transform
through :class:`preprocessing.PowerTransformer` to map data
from various distributions to a Gaussian distribution.

Box-Cox is useful as a transformation in modeling problems
where homoscedasticity and normality are desired. Below
are examples of Box-Cox applied to three distributions: a
heavily left-skewed lognormal distribution, a left-skewed
chi-squared distribution, and a right-skewed Weibull
distribution.

Note that when Box-Cox is applied to the lognormal distribution,
the estimated parameter (lambda) takes on a value close to
zero, since Box-Cox is defined as the natural logarithm when
lambda = 0. In all cases, the maximum likelihood estimated
lambdas achieve very close results to the normal distribution.
"""
print(__doc__)

# Author: Eric Chang <ericchang2017@u.northwestern.edu>
# License: BSD 3 clause

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from sklearn.preprocessing import PowerTransformer

N_SAMPLES = 3000
SEED = 304

pt = PowerTransformer(method='box-cox')
rng = np.random.RandomState(SEED)
size = (N_SAMPLES, 1)


# lognormal distribution
X_lognormal = rng.lognormal(size=size)

# chi-squared distribution
df = 3
X_chisq = rng.chisquare(df=df, size=size)

# weibull distribution
a = 50
X_weibull = rng.weibull(a=a, size=size)

# gaussian distribution
loc = 100
X_gaussian = rng.normal(loc=loc, size=size)

# uniform distirbution
X_uniform = rng.uniform(low=0, high=1, size=size)

# bimodal distribution
loc_a, loc_b = 100, 105
X_a, X_b = rng.normal(loc=loc_a, size=size), rng.normal(loc=loc_b, size=size)
X_bimodal = np.concatenate([X_a, X_b], axis=0)


# create plots
distributions = [
    ('Lognormal', X_lognormal),
    ('Chi-squared', X_chisq),
    ('Weibull', X_weibull),
    ('Gaussian', X_gaussian),
    ('Uniform', X_uniform),
    ('Bimodal', X_bimodal),
]

colors = ['firebrick', 'darkorange', 'goldenrod',
          'seagreen', 'royalblue', 'darkorchid']

fig, axes = plt.subplots(nrows=4, ncols=3)
axes_idxs = [(0, 3), (1, 4), (2, 5), (6, 9), (7, 10), (8, 11)]
axes = axes.flatten().tolist()

params = {'font.size': 6, 'hist.bins': 150}
matplotlib.rcParams.update(params)

for distribution, color, idxs in zip(distributions, colors, axes_idxs):
    name, X = distribution
    ax_original, ax_trans = axes[idxs[0]], axes[idxs[1]]

    # perform power transform
    X_trans = pt.fit_transform(X)
    lmbda = round(pt.lambdas_[0], 2)

    ax_original.hist(X, color=color)
    ax_original.set_title(name)
    ax_trans.hist(X_trans, color=color)
    ax_trans.set_title('{} after Box-Cox, $\lambda$ = {}'.format(name, lmbda))

plt.tight_layout()
plt.show()
