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

X_min_max_scaled = MinMaxScaler().fit_transform(X)
X_max_abs_scaled = MaxAbsScaler().fit_transform(X)
X_standard_scaled = StandardScaler().fit_transform(X)
X_robust_scaled = RobustScaler(quantile_range=(25, 75)).fit_transform(X)
X_l2_normalized = Normalizer().fit_transform(X)
X_quantile_normalized = QuantileNormalizer().fit_transform(X)

y = minmax_scale(y_full)  # To make colors corresponding to the target


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

fig = plt.figure(figsize=(15, 50))
gs = gridspec.GridSpec(15, 5,
                       width_ratios=[5, 1, 0.1, 5, 1], wspace=0.3,
                       height_ratios=[5, 1] * 7 + [0.4], hspace=0.4)
subplots = list(plt.subplot(g) for g in gs)

for i, (X, title) in enumerate((
        (X, "Unscaled data"),
        (X_min_max_scaled, "Data after min-max scaling"),
        (X_robust_scaled, "Data after robust scaling"),
        (X_max_abs_scaled, "Data after max-abs scaling"),
        (X_standard_scaled, "Data after standard scaling"),
        (X_l2_normalized, "Data after sample-wise L2 normalizing"),
        (X_quantile_normalized, "Data after quantile normalizing"))):
    offset = 10 * i

    # Distribution with all outliers
    axes = subplots[offset:offset + 2] + subplots[offset + 5:offset + 7]
    plot_distribution(axes, X, y, hist_nbins=50,
                      plot_title=title + " including outliers\n",
                      X_label="Median Income", y_label="Number of households")

    # Some blank vertical space between two plots so they don't overlap
    subplots[offset + 2].axis('off')
    subplots[offset + 7].axis('off')

    # Distribution with extreme outliers removed
    X0_min, X0_99th_pc = np.percentile(X[:, 0], [0, 99])
    X1_min, X1_99th_pc = np.percentile(X[:, 1], [0, 99])

    non_outliers = np.all(X < [X0_99th_pc, X1_99th_pc], axis=1)
    axes = subplots[offset + 3:offset + 5] + subplots[offset + 8:offset + 10]
    plot_distribution(axes, X[non_outliers], y[non_outliers], hist_nbins=50,
                      plot_title=(title +
                                  "\nZoomed-in at percentile range [0, 99)"),
                      X_label="Median Income", y_label="Number of households")

# Plot a heatmap legend for the y, combining a row of 4 cols
heatmap_legend_ax = plt.subplot(gs[-5:])
norm = mpl.colors.Normalize(y_full.min(), y_full.max())
mpl.colorbar.ColorbarBase(heatmap_legend_ax, cmap=cm.plasma_r,
                          norm=norm, orientation='horizontal',
                          label='Color mapping for values of y')
plt.show()
