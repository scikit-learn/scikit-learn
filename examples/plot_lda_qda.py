"""
==============================================================
Linear Discriminant Analysis & Quadratic Discriminant Analysis
==============================================================

Plot the confidence ellipsoids of each class and decision boundary
"""

import numpy as np
import pylab as pl
import matplotlib as mpl
from matplotlib import collections

from scikits.learn.lda import LDA

################################################################################
# generate datasets
def dataset_fixed_cov():
    '''Generate 2 Gaussians samples with the same covariance matrix'''
    n, dim = 300, 2
    np.random.seed(0)
    C = np.array([[0., -0.7], [3.5, .7]])
    X = np.r_[np.dot(np.random.randn(n, dim), C),
              np.dot(np.random.randn(n, dim), C) + np.array([3, 3])]
    y = np.hstack((np.zeros(n), np.ones(n)))
    return X, y

def dataset_cov():
    '''Generate 2 Gaussians samples with different covariance matrices'''
    n, dim = 300, 2
    np.random.seed(0)
    C = np.array([[0., -0.7], [3.5, .7]])
    X = np.r_[np.dot(np.random.randn(n, dim), C),
              np.dot(np.random.randn(n, dim), C.T) + np.array([3, 3])]
    y = np.hstack((np.zeros(n), np.ones(n)))
    return X, y

################################################################################
# boundaries fonction

def lda_decision_boundary(lda, x):
    '''
    lda : lda classifier instance.
    x :   abscissa coordinate.

    Return ordinate coordinate.
    '''
    m0, m1 = lda.means_
    p0, p1 = lda.priors_
    S = lda.covariance_
    Sinv = np.linalg.inv(S)
    z = np.dot(Sinv, m0 - m1)
    a = 0.5 * np.dot((m0 + m1).T, z) - np.log(p0 / p1)
    return (a - x * z[0]) / z[1]

################################################################################
# plot functions

def lda_plot(lda, X, y, y_pred, fig_index):
    splot = pl.subplot(2, 2, fig_index)
    if fig_index == 1:
        pl.title('Linear Discriminant Analysis')
        pl.ylabel('Fixed covariance')
    elif fig_index == 3:
        pl.ylabel('Different covariances')

    tp = (y == y_pred) # True Positive
    tp0, tp1 = tp[y == 0], tp[y == 1]
    X0, X1 = X[y == 0], X[y == 1]
    X0_tp, X0_fp = X0[tp0], X0[tp0 != True]
    X1_tp, X1_fp = X1[tp1], X1[tp1 != True]
    xmin, xmax = X[:, 0].min(), X[:, 0].max()
    ymin, ymax = X[:, 1].min(), X[:, 1].max()
    b_xmin, b_xmax = lda_decision_boundary(lda, np.array([xmin, xmax]))

    # class 0: area
    pl.fill([xmin, xmax, xmax, xmin, xmin],
            [b_xmin, b_xmax, ymin, ymin, b_xmin],
            facecolor='#ffbbbb', edgecolor='none') # light red

    # class 1: area
    pl.fill([xmin, xmax, xmax, xmin, xmin],
            [b_xmin, b_xmax, ymax, ymax, b_xmin],
            facecolor='#bbbbff', edgecolor='none') # light blue

    # class 0: dots 
    pl.plot(X0_tp[:, 0], X0_tp[:, 1], 'o', color='red')
    pl.plot(X0_fp[:, 0], X0_fp[:, 1], '.', color='#990000') # dark red

    # class 1: dots
    pl.plot(X1_tp[:, 0], X1_tp[:, 1], 'o', color='blue')
    pl.plot(X1_fp[:, 0], X1_fp[:, 1], '.', color='#000099') # dark blue

    # means
    pl.plot(lda.means_[0][0], lda.means_[0][1],
            'o', color = 'black', markersize=10)
    pl.plot(lda.means_[1][0], lda.means_[1][1],
            'o', color = 'black', markersize=10)

    # decision boundary
    pl.plot([xmin, xmax], [b_xmin, b_xmax], '-', color='black')

    # covariances
    v, w = np.linalg.eigh(lda.covariance_)
    u = w[0] / np.linalg.norm(w[0])
    angle = np.arctan(u[1]/u[0])
    angle = 180 * angle / np.pi # convert to degrees
    ell = mpl.patches.Ellipse(lda.means_[0], v[0], v[1],
                                180 + angle, color='red')
    ell.set_clip_box(splot.bbox)
    ell.set_alpha(0.5)
    splot.add_artist(ell)
    ell = mpl.patches.Ellipse (lda.means_[1], v[0], v[1],
                                180 + angle, color='blue')
    ell.set_clip_box(splot.bbox)
    ell.set_alpha(0.5)
    splot.add_artist(ell)

    pl.axis('tight')

def main():
    for i, (X, y) in enumerate([dataset_fixed_cov(), dataset_cov()]):
        # LDA
        lda = LDA()
        y_pred = lda.fit(X, y, store_covariance=True).predict(X)
        lda_plot(lda, X, y, y_pred, 2 * i + 1)

        # QDA #FIXME
        lda = LDA()
        y_pred = lda.fit(X, y, store_covariance=True).predict(X)
        lda_plot(lda, X, y, y_pred, 2 * i + 2)
    pl.suptitle('LDA vs QDA')
    pl.show()

if __name__ == '__main__': main()
