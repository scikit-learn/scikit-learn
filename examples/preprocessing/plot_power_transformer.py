"""
==========================================================
Using PowerTransformer to apply the Box-Cox transformation
==========================================================

This example demonstrates the use of the Box-Cox transform through
:class:`preprocessing.PowerTransformer` to map data from various distributions
to a normal distribution.

Box-Cox is useful as a transformation in modeling problems where
homoscedasticity and normality are desired. Below are examples of Box-Cox
applied to six different probability distributions: Lognormal, Chi-squared,
Weibull, Gaussian, Uniform, and Bimodal.

Note that the transformation successfully maps the data to a normal
distribution when applied to certain datasets, but is ineffective with others.
This highlights the importance of visualizing the data before and after
transformation. Also note that while the standardize option is set to False for
the plot examples, by default, :class:`preprocessing.PowerTransformer` also
applies zero-mean, unit-variance standardization to the transformed outputs.
"""

# Author: Eric Chang <ericchang2017@u.northwestern.edu>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import PowerTransformer, minmax_scale

print(__doc__)


N_SAMPLES = 3000
FONT_SIZE = 6
BINS = 100


pt = PowerTransformer(method='box-cox', standardize=False)
rng = np.random.RandomState(304)
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

# uniform distribution
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
    ('Bimodal', X_bimodal)
]

colors = ['firebrick', 'darkorange', 'goldenrod',
          'seagreen', 'royalblue', 'darkorchid']

fig, axes = plt.subplots(nrows=4, ncols=3)
axes = axes.flatten()
axes_idxs = [(0, 3), (1, 4), (2, 5), (6, 9), (7, 10), (8, 11)]
axes_list = [(axes[i], axes[j]) for i, j in axes_idxs]


for distribution, color, axes in zip(distributions, colors, axes_list):
    name, X = distribution
    # scale all distributions to the range [0, 10]
    X = minmax_scale(X, feature_range=(1e-10, 10))

    # perform power transform
    X_trans = pt.fit_transform(X)
    lmbda = round(pt.lambdas_[0], 2)

    ax_original, ax_trans = axes

    ax_original.hist(X, color=color, bins=BINS)
    ax_original.set_title(name, fontsize=FONT_SIZE)
    ax_original.tick_params(axis='both', which='major', labelsize=FONT_SIZE)

    ax_trans.hist(X_trans, color=color, bins=BINS)
    ax_trans.set_title('{} after Box-Cox, $\lambda$ = {}'.format(name, lmbda),
                       fontsize=FONT_SIZE)
    ax_trans.tick_params(axis='both', which='major', labelsize=FONT_SIZE)


plt.tight_layout()
plt.show()
