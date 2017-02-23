#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
=============================================================
Compare the effect of different scalers on data with outliers
=============================================================

The feature 0 and feature 5 of california housing dataset contains large
outliers that can make visualization of the data difficult.

Also linear models like :class:`sklearn.linear_model.SVM` require data which is
approximately normalized to the [-1, 1] or [0, 1] range, or at the very least
have all the features on the same scale.

This example uses different scalers and normalizers to bring the data within a
smaller range.
"""
from __future__ import print_function
print(__doc__)

# Author:  Raghav RV <rvraghav93@gmail.com>
#          Thomas Unterthiner
# License: BSD 3 clause

from collections import OrderedDict

import numpy as np

import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import cm, gridspec

from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import minmax_scale
from sklearn.preprocessing import MaxAbsScaler
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import Normalizer
from sklearn.preprocessing.data import QuantileNormalizer

from sklearn.datasets import fetch_california_housing

dataset = fetch_california_housing()
X_full, y_full = dataset.data, dataset.target

# Take only 2 features to make visualization easier
# Feature of 0 has a long tail distribution.
# Feature 5 has a few but very large outliers.

X = X_full[:, [0, 5]]

distributions = OrderedDict((
    ('Unscaled data', X),
    ('Data after min-max scaling',
        MinMaxScaler().fit_transform(X)),
    ('Data after robust scaling',
        RobustScaler(quantile_range=(25, 75)).fit_transform(X)),
    ('Data after max-abs scaling',
        MaxAbsScaler().fit_transform(X)),
    ('Data after standard scaling',
        StandardScaler().fit_transform(X)),
    ('Data after sample-wise L2 normalizing',
        Normalizer().fit_transform(X)),
    ('Data after quantile normalizing (uniform pdf)',
        QuantileNormalizer(output_pdf='uniform')
        .fit_transform(X)),
    ('Data after quantile normalizing (gaussian pdf)',
        QuantileNormalizer(output_pdf='norm')
        .fit_transform(X))))

y = minmax_scale(y_full)  # To make colors corresponding to the target),


def plot_distribution(axes, X, y, hist_nbins=50, plot_title="", size=(15, 10),
                      X_label="", y_label=""):
    ax, hist_X1, hist_X0, empty = axes
    empty.axis('off')

    ax.set_title(plot_title, fontsize=12)
    ax.set_xlabel(X_label)
    ax.set_ylabel(y_label)

    # The scatter plot
    colors = cm.plasma_r(y)
    ax.scatter(X[:, 0], X[:, 1], alpha=0.5, marker='o', s=5, lw=0, c=colors)

    # Removing the top and the right spine for aesthetics
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # Histogram for axis X1 (feature 5)
    hist_X1.set_ylim(ax.get_ylim())
    hist_X1.hist(X[:, 1], bins=hist_nbins, orientation='horizontal',
                 color='grey', ec='grey')
    hist_X1.axis('off')

    # Histogram for axis X0 (feature 0)
    hist_X0.set_xlim(ax.get_xlim())
    hist_X0.invert_yaxis()
    hist_X0.hist(X[:, 0], bins=hist_nbins, orientation='vertical',
                 color='grey', ec='grey')
    hist_X0.axis('off')

n_dist = len(distributions)
fig = plt.figure(figsize=(15, n_dist * 8 + 1))
gs = gridspec.GridSpec(n_dist * 2 + 1, 5,
                       width_ratios=[5, 1, 0.1, 5, 1], wspace=0.3,
                       height_ratios=[5, 1] * n_dist + [0.4],
                       hspace=0.4)
subplots = list(plt.subplot(g) for g in gs)

for i, (title, X) in enumerate(distributions.items()):
    offset = 10 * i
    # Distribution with all outliers
    axes = subplots[offset:offset + 2] + subplots[offset + 5:offset + 7]
    plot_distribution(axes, X, y, hist_nbins=200,
                      plot_title=title + " including outliers\n",
                      X_label="Median Income", y_label="Number of households")

    # Some blank vertical space between two plots so they don't overlap
    subplots[offset + 2].axis('off')
    subplots[offset + 7].axis('off')

    zoom_in_percentile_range = (0, 99)
    # Distribution with extreme outliers removed
    cutoffs_X0 = np.percentile(X[:, 0], zoom_in_percentile_range)
    cutoffs_X1 = np.percentile(X[:, 1], zoom_in_percentile_range)

    non_outliers_mask = (
        np.all(X > [cutoffs_X0[0], cutoffs_X1[0]], axis=1) &
        np.all(X < [cutoffs_X0[1], cutoffs_X1[1]], axis=1))
    axes = subplots[offset + 3:offset + 5] + subplots[offset + 8:offset + 10]
    plot_distribution(axes, X[non_outliers_mask], y[non_outliers_mask],
                      hist_nbins=50,
                      plot_title=(title +
                                  "\nZoomed-in at percentile range %s"
                                  % str(zoom_in_percentile_range)),
                      X_label="Median Income", y_label="Number of households")

# Plot a heatmap legend for the y, combining a row of 4 cols
heatmap_legend_ax = plt.subplot(gs[-5:])
norm = mpl.colors.Normalize(y_full.min(), y_full.max())
mpl.colorbar.ColorbarBase(heatmap_legend_ax, cmap=cm.plasma_r,
                          norm=norm, orientation='horizontal',
                          label='Color mapping for values of y')
plt.show()
