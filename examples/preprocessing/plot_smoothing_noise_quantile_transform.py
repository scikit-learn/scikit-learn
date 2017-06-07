#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
========================================================
Effect of smoothing noise when using QuantileTransformer
========================================================

The parameter ``smoothing_noise`` can be used if some specific feature values
are repeated identically many times to the point of being predominant in the
dataset. This is typically be observed when the feature encode ordinal
values such as user ratings, prices, coarse-grained units of time, etc. By
default, a small Gaussian noise is added during ``fit`` time.

Without smoothing noise, the ``QuantileTransformer`` will map those values to
some arbitrary value: the highest quantile value for all the inputs with the
same value. While this is usually not an issue when ``QuantileTransformer`` is
used as a preprocessing transformer for a subsequent supervised estimator, it
can lead to surprising results when manually inspecting the transformed values
(e.g. for visualization or reporting).

The goal of the ``smoothing_noise`` is to make it possible to map those
repeated values to some middle quantile value to make interpretation more
intuitive as demonstrated in the following.

"""

# Author:  Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import QuantileTransformer

print(__doc__)

N_QUANTILES = 1000
FEAT_VAL = 3.0


def plot_transform_feat_val(ax, transformer, title):
    """Plot the mapping function as well as a specific feature value."""
    ref = np.linspace(0, 1, num=N_QUANTILES)

    ax.plot(transformer.quantiles_, ref)
    ax.scatter(FEAT_VAL, transformer.transform(FEAT_VAL), c='r',
               label=r'$f({0}) = {1:.2f}$'.format(
                   FEAT_VAL,
                   np.ravel(transformer.transform(FEAT_VAL))[0]))
    ax.set_xlabel('Feature values')
    ax.set_ylabel('Quantiles in %')
    ax.set_title(title)
    ax.legend(loc=4)
    # make nice axis layout
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlim([1, 5.1])
    ax.set_ylim([0, 1])
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))


###############################################################################
# We can create a synthetic dataset representing the customers'
# ratings for a restaurant. The scale used is ranging from 1 to 5 and
# a large number of customers attributed a grade of 3 to the current
# restaurant.
#
# The ``smoothing_noise`` can be disabled in ``QuantileTransformer``.
# When dealing with a data set with a predominant value, this feature
# value can be affected to several quantiles. When provided to the transformer,
# this feature value will be mapped to the largest quantile. In practice,
# machine learning algorithms will usually not be affected by such
# characteristics. However, manual interpretation might be counter intuitive.
#
# From the below plot, we would expect that a vote corresponding to
# the value 3 would be mapped to the median (e.g., 0.5). However, the
# default behaviour of the 'interp' numpy function will map this
# feature value to the greater quantile as shown by the marker in the
# figure.
#
# A solution is to apply a small smoothing noise before computing the
# quantiles. The parameter ``smoothing_noise=True`` (default behaviour) offers
# this possibility as illustrated above.  In this case, the marker is centered
# at the median as expected.

X = np.array([1] * 2000 +
             [2] * 1000 +
             [3] * 7000 +
             [4] * 2000 +
             [5] * 1000).reshape(-1, 1)

# create the subplots
_, (ax1, ax2) = plt.subplots(1, 2)

qt = QuantileTransformer(n_quantiles=N_QUANTILES,
                         smoothing_noise=False)
qt.fit(X)
plot_transform_feat_val(ax1, qt, 'Without smoothing')

qt = QuantileTransformer(n_quantiles=N_QUANTILES,
                         smoothing_noise=True)
qt.fit(X)
plot_transform_feat_val(ax2, qt, 'With smoothing')
plt.tight_layout()
plt.show()
