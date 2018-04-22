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
transformation. We set the standardize option is set to True for these examples,
implying that we have performed a zero-mean, unit-variance standardization to
the transformed outputs.

We also include overlays of transformations via
:class:`preprocessing.QuantileTransformer` for comparison purposes. The blue
scatterplots represent the quantile transformation, while the red scatterplots
represent the power transformation.
"""

# Authors: Theodore Yoong <guo-zheng.yoong@univ.ox.ac.uk>, Eric Chang <ericchang2017@u.northwestern.edu>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import QuantileTransformer, PowerTransformer

print(__doc__)

qt = QuantileTransformer(output_distribution='normal')
pt = PowerTransformer(method='box-cox', standardize=True)


# lognormal distribution
X_lognormal = np.random.lognormal(0, 1, 1000)

# chi-squared distribution
X_chisq = np.random.chisquare(2,1000)

# weibull distribution
X_weibull = np.random.weibull(5, 1000)

# gaussian distribution
X_gaussian = np.random.normal(0,1,1000)
# ensure strictly positive values (box-cox does not allow non-positive values)
X_gaussian = X_gaussian - X_gaussian.min() + 1

# uniform distribution
X_uniform = np.random.uniform(0,1,1000)

# bimodal distribution
X_bimodal = 0.5*(np.random.normal(0,1,1000) + np.random.normal(0,50,1000))
# ensure strictly positive values
X_bimodal = X_bimodal - X_bimodal.min() + 1


# create plots
distributions = [
    ('Lognormal', X_lognormal),
    ('Chi-squared', X_chisq),
    ('Weibull', X_weibull),
    ('Gaussian', X_gaussian),
    ('Uniform', X_uniform),
    ('Bimodal', X_bimodal)
]

fig = plt.figure(figsize = (10,10))

for distribution, i in zip(distributions, range(1,7)):
    name, X = distribution
    # transform X (currently a row vector) into a column vector
    # before performing respective transforms
    X_qt_trans = qt.fit_transform(X.reshape(-1,1))
    X_pt_trans = pt.fit_transform(X.reshape(-1,1))
    
    ax = fig.add_subplot(3,2,i)
    ax.set_title('{}'.format(name))
    # we denote 's' as markersize, since scatter does not have markersize kwarg
    plt.scatter(X, X_qt_trans, s=1)
    plt.scatter(X, X_pt_trans, s=1, c='r')

plt.tight_layout()
plt.show()
