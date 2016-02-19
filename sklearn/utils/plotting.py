"""
Utility function for plotting the decision regions of a classifier.
"""

# Authors: Sebastian Raschka <mail@sebastianraschka.com.com>
#
# Licence: BSD 3 clause


# Sebastian Raschka 2014-2016
# mlxtend Machine Learning Library Extensions
#
# A function for plotting decision regions of classifiers.
# Author: Sebastian Raschka <sebastianraschka.com>
#
# License: BSD 3 clause


from itertools import cycle
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap


def plot_decision_regions(X, y, clf,
                          ax=None,
                          X_highlight=None,
                          res=0.02, legend=1,
                          hide_spines=True,
                          markers='s^oxv<>',
                          colors='red,blue,limegreen,gray,cyan'):
    """Plot decision regions of a classifier.

    Parameters
    ----------
    X : array-like, shape = [n_samples, n_features]
        Feature Matrix.
    y : array-like, shape = [n_samples]
        True class labels.
    clf : Classifier object.
        Must have a .predict method.
    ax : matplotlib.axes.Axes (default: None)
        An existing matplotlib Axes. Creates
        one if ax=None.
    X_highlight : array-like, shape = [n_samples, n_features] (default: None)
        An array with data points that are used to highlight samples in `X`.
    res : float (default: 0.02)
        Grid width. Lower values increase the resolution but
        slow down the plotting.
    hide_spines : bool (default: True)
        Hide axis spines if True.
    legend : int (default: 1)
        Integer to specify the legend location.
        No legend if legend is 0.
    markers : list
        Scatterplot markers.
    colors : str (default 'red,blue,limegreen,gray,cyan')
        Comma separated list of colors.

    Returns
    ---------
    ax : matplotlib.axes.Axes object

    """
    # check if data is numpy array
    for a in (X, y):
        if not isinstance(a, np.ndarray):
            raise ValueError('%s must be a NumPy array.' % a.__name__)

    if ax is None:
        ax = plt.gca()

    # check if test data is provided
    plot_testdata = True
    if not isinstance(X_highlight, np.ndarray):
        if X_highlight is not None:
            raise ValueError('X_test must be a NumPy array or None')
        else:
            plot_testdata = False

    if len(X.shape) == 2 and X.shape[1] > 1:
        dim = '2d'
    else:
        dim = '1d'

    marker_gen = cycle(list(markers))

    # make color map
    n_classes = len(np.unique(y.astype(int)))
    colors = colors.split(',')
    cmap = ListedColormap(colors[:n_classes])

    # plot the decision surface
    if dim == '2d':
        y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    else:
        y_min, y_max = -1, 1

    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, res),
                         np.arange(y_min, y_max, res))

    if dim == '2d':
        y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
        Z = clf.predict(np.array([xx.ravel(), yy.ravel()]).T)
    else:
        y_min, y_max = -1, 1
        Z = clf.predict(np.array([xx.ravel()]).T)

    Z = Z.reshape(xx.shape)
    ax.contourf(xx, yy, Z, alpha=0.3, cmap=cmap)

    ax.axis(xmin=xx.min(), xmax=xx.max(), y_min=yy.min(), y_max=yy.max())

    # plot class samples

    for c in np.unique(y.astype(int)):
        if dim == '2d':
            y_data = X[y.astype(int) == c, 1]
        else:
            y_data = [0 for i in X[y.astype(int) == c]]

        ax.scatter(x=X[y.astype(int) == c, 0],
                   y=y_data,
                   alpha=0.8,
                   c=cmap(c),
                   marker=next(marker_gen),
                   label=c)

    if hide_spines:
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    if not dim == '2d':
        ax.axes.get_yaxis().set_ticks([])

    if legend:
        legend = plt.legend(loc=legend,
                            fancybox=True,
                            framealpha=0.3,
                            scatterpoints=1,
                            borderaxespad=1)

        ax.add_artist(legend)

    if plot_testdata:
        if dim == '2d':
            ax.scatter(X_highlight[:, 0],
                       X_highlight[:, 1],
                       c='',
                       alpha=1.0,
                       linewidth=1,
                       marker='o',
                       s=80)
        else:
            ax.scatter(X_highlight,
                       [0 for i in X_highlight],
                       c='',
                       alpha=1.0,
                       linewidth=1,
                       marker='o',
                       s=80)

    return ax
