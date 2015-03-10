# Authors: Cecilia Damon <cecilia.damon@institut-hypercube.org>

import numpy as np
import matplotlib.pylab as MP

from sklearn import datasets
from sklearn.utils.testing import assert_array_equal
from sklearn.feature_selection import (chi2, f_classif, f_regression,GenericUnivariateSelect)
from sklearn.feature_selection.MI_RenyiCC import (_MI_Filter,univariate_f_MI,univariate_forward_f_MI,\
                                                  multivariate_backward_f_MI,multivariate_forward_f_MI, MI_RenyiCC_Multi)

def test_equal_dependency_bivariate_data():
    '''
    Three examples of bivariate variable (x,y) testing both their perfect dependency (with x=y) vs
    independency (random generation of x and y). Each example considers the different possible
    data type configurations, i.e. x and y are discrete ('dd'), x and y are continuous ('cc') and
    x is continuous while y is discrete ('cd')
    '''
    x= np.random.randint(0,3,10)
    y= np.random.randint(0,10, 10)
    y2 = x
    print("MI score of two discrete random variables:", MI_RenyiCC_Multi(np.vstack((x,y)).T, type="d"))
    print("MI score of two discrete dependent variables:", MI_RenyiCC_Multi(np.vstack((x,y2)).T, type="d"))
    MP.figure()
    MP.plot(x,y,'r.')
    MP.plot(x,y2,'b.')

    x= np.random.randn(100)
    y= np.random.randn(100)
    y2 = x
    print("MI score of two continuous random variables:", MI_RenyiCC_Multi(np.vstack((x,y)).T, type="c"))
    print("MI score of two continuous dependent variables:", MI_RenyiCC_Multi(np.vstack((x,y2)).T, type="c"))
    MP.figure()
    MP.plot(x,y,'r.')
    MP.plot(x,y2,'b.')

    x= 5*np.random.randn(100)+5
    y= np.random.randint(0,10, 100)
    y2 = np.asarray(np.round(x),int)
    print("MI score of a pair of (continuous, discrete) random variables:", MI_RenyiCC_Multi(x.reshape((100,1)),y, type="cd"))
    print("MI score of a pair of (continuous, discrete) dependent variables:", MI_RenyiCC_Multi(x.reshape((100,1)),y2, type="cd"))
    MP.figure()
    MP.plot(x,y,'r.')
    MP.plot(x,y2,'b.')

def rotate(phi, x):
    """
    Rotate 2D points by angle phi
    """
    from math import cos, sin
    from numpy import inner
    z = [[cos(phi), sin(phi)], [-sin(phi), cos(phi)]]
    return inner(x, z)

def test_correlation_examples(N=500):
    '''
    Python code of correlation examples from beaucronin : https://gist.github.com/beaucronin/2509755
    Python translation of examples http://en.wikipedia.org/wiki/File:Correlation_examples2.svg
    Title: An example of the correlation of x and y for various distributions of (x,y) pairs
    Author: Denis Boigelot

    Parameters
    ----------
    N : the number of samples

    Returns
    -------

    '''

    from numpy.random import (
    uniform as runif,
    multivariate_normal as rmvn,
    normal as rnorm
    )
    from numpy import inner, linspace, array, vstack
    from math import pi, cos, sin, pow, sqrt
    import matplotlib.pyplot as plt

    datasets = []
    MI_scores = []
    for corr in [1., .8, .4, 0., -.4, -.8, -1.]:
        x = rmvn([0., 0.], [[1., corr], [corr, 1.]], N)
        mi = MI_RenyiCC_Multi(x, type="c")
        print("MI score of two continuous linear correlated (correlation degree of %f) gaussian variables : %f" % (corr, mi))
        MI_scores.append(mi)
        datasets.append(x)

    for phi in [0., pi/12., pi/6., pi/4., pi/2. - pi/6., pi/2. - pi/12, pi/2]:
        x = rmvn([0., 0.], [[1., 1.], [1., 1.]], N)
        x = rotate(phi, x)
        mi = MI_RenyiCC_Multi(x, type="c")
        print("MI score of two continuous linear correlated (correlation slope of %f) gaussian variables : %f" % (phi, mi))
        MI_scores.append(mi)
        datasets.append(x)

    a = linspace(-1, 1, N)
    x = array([(x0, 4. * pow(x0 * x0 - .5, 2.) + runif(-1./3., 1./3., 1))
        for x0 in a])
    mi = MI_RenyiCC_Multi(x, type="c")
    print("MI score of two continuous non-linear correlated (y=4*(x0^2-0.5)^2+c) variables : %f" %  mi)
    MI_scores.append(mi)
    datasets.append(x)

    x = rotate(-pi/8., array([(x0, runif(-1., 1.)) for x0 in a]))
    mi = MI_RenyiCC_Multi(x, type="c")
    print("MI score of two continuous non-linear correlated variables : %f" %  mi)
    MI_scores.append(mi)
    datasets.append(x)

    x = rotate(-pi/8, x)
    mi = MI_RenyiCC_Multi(x, type="c")
    print("MI score of two continuous non-linear correlated variables : %f" %  mi)
    MI_scores.append(mi)
    datasets.append(x)

    x = array([(x0, x0 * x0 + runif(-.5, .5)) for x0 in a])
    mi = MI_RenyiCC_Multi(x, type="c")
    print("MI score of two continuous non-linear correlated variables : %f" %  mi)
    MI_scores.append(mi)
    datasets.append(x)

    signs = [1. if runif() < .5 else -1. for _ in range(N)]
    x = array([(x0, (x0 * x0 + runif(0., .5)) * sign)
        for x0, sign in zip(a, signs)])
    mi = MI_RenyiCC_Multi(x, type="c")
    print("MI score of two continuous non-linear correlated variables : %f" %  mi)
    MI_scores.append(mi)
    datasets.append(x)

    x = array([(sin(x0 * pi) + rnorm(0., .125), cos(x0 * pi) + rnorm(0., .125))
        for x0 in a])
    mi = MI_RenyiCC_Multi(x, type="c")
    print("MI score of two continuous non-linear correlated variables : %f" %  mi)
    MI_scores.append(mi)
    datasets.append(x)

    x = vstack((
        rmvn([3., 3], [[1., 0.], [0., 1.]], round(N/4)),
        rmvn([-3., 3], [[1., 0.], [0., 1.]], round(N/4)),
        rmvn([-3., -3], [[1., 0.], [0., 1.]], round(N/4)),
        rmvn([3., -3], [[1., 0.], [0., 1.]], round(N/4))
        ))
    mi = MI_RenyiCC_Multi(x, type="c")
    print("MI score of two continuous non-linear correlated variables : %f" %  mi)
    MI_scores.append(mi)
    datasets.append(x)

    """
    Plot the datasets, mimicking the original plot from Wikipedia.
    """
    plt.figure()
    print("MI_scores:",MI_scores)
    for i in range(len(datasets)):
        plt.subplot(3, 7, i+1)
        x = [a[0] for a in datasets[i]]
        y = [a[1] for a in datasets[i]]
        plt.plot(x, y, '.', markersize=1.)
        plt.title("%.4f" %MI_scores[i])
        # plt.axis('scaled')
        plt.xticks([])
        plt.yticks([])
        ax = plt.gca()
        ax.set_axis_off()
        if i == 14:
            plt.xlim([-1, 1])
            plt.ylim([-1./3., 1.+1./3.])
        elif i == 15:
            z = sqrt(2. + sqrt(2.)) / sqrt(2.)
            plt.xlim([-z, z])
            plt.ylim([-z, z])
        elif i == 16:
            plt.xlim([-sqrt(2.), sqrt(2.)])
            plt.ylim([-sqrt(2.), sqrt(2.)])
        elif i == 17:
            plt.xlim([-1, 1])
            plt.ylim([-.5, 1.5])
        elif i == 18:
            plt.xlim([-1.5, 1.5])
            plt.ylim([-1.5, 1.5])
        elif i == 19:
            plt.xlim([-1.5, 1.5])
            plt.ylim([-1.5, 1.5])
        elif i == 20:
            plt.xlim([-7, 7])
            plt.ylim([-7, 7])
        else:
            plt.xlim([-4, 4])
            plt.ylim([-4, 4])
        ax.set_aspect('equal', adjustable='datalim')
    #plt.savefig('out.pdf')

def test_skdatasets_classif_MIFilter_vs_f_classif(X, y, p, k):
    '''
    Test several scikits learn datasets to compare features selection based on MI and anova
    for classification problems
    '''
    #scores, r = univariate_f_MI(X, y, type='cd',njobs=4)
    #scores, r= univariate_f_MI(X, y,k=30, type='cd') njobs=4
    #print("MI scores:",scores)
    #support_MI = get_support(scores, 4)
    kbest_univ_f_MI = _MI_Filter(univariate_f_MI, mode='k_best',param=k,type='cd',njobs=4)
    kbest_univ_f_MI.fit(X,y)
    support_MI = kbest_univ_f_MI._get_support_mask()
    print("MI Univariate scores:",kbest_univ_f_MI.scores_)
    print("MI Univariate r:",kbest_univ_f_MI.ranking_)
    print("MI Univariate support:",support_MI)
    kfirst_uforward_MI = _MI_Filter(univariate_forward_f_MI, mode='k_first',param=k,type='cd',njobs=4)
    kfirst_uforward_MI.fit(X,y)
    support_uforward_MI = kfirst_uforward_MI._get_support_mask()
    print("MI Univariate forward scores:",kfirst_uforward_MI.scores_)
    print("MI Univariate forward r:",kfirst_uforward_MI.ranking_)
    print("MI Univariate forward support:",support_uforward_MI)
    kfirst_mforward_MI = _MI_Filter(multivariate_forward_f_MI, mode='k_first',param=k,type='cd',njobs=4)
    kfirst_mforward_MI.fit(X,y)
    support_mforward_MI = kfirst_mforward_MI._get_support_mask()
    print("MI multivariate forward scores:",kfirst_mforward_MI.scores_)
    print("MI multivariate forward r:",kfirst_mforward_MI.ranking_)
    print("MI multivariate forward support:",support_mforward_MI)

    kfirst_mbackward_MI = _MI_Filter(multivariate_backward_f_MI, mode='k_first',param=-k,type='cd',njobs=4)
    kfirst_mbackward_MI.fit(X,y)
    support_mbackward_MI = kfirst_mbackward_MI._get_support_mask()
    print("MI multivariate backward scores:",kfirst_mbackward_MI.scores_)
    print("MI multivariate backward r:",kfirst_mbackward_MI.ranking_)
    print("MI multivariate backward support:",support_mbackward_MI)
    filter_F = GenericUnivariateSelect(f_classif, mode='k_best',param=k)
    filter_F.fit(X, y)
    print("F scores:",filter_F.scores_)
    support_F = filter_F._get_support_mask()
    print("F support :",support_F)
    '''
    gtruth = np.zeros(p+k)
    gtruth[:k] = 1
    print("gtruth:",gtruth)
    assert_array_equal(support_MI, gtruth)
    assert_array_equal(support_F, gtruth)
    assert_array_equal(support_uforward_MI, gtruth)
    assert_array_equal(support_mforward_MI, gtruth)
    assert_array_equal(support_mbackward_MI, gtruth)
    '''

def test_iris_classif_MIFilter_vs_f_classif():
    '''y is a Linear combination of 4 features among 14. Compare informative features obtained with MI to those
    obtained with Ftest'''
    print("Iris Dataset\n")
    #p=16; k =4
    p=3; k=4
    data = datasets.load_iris()
    X = data.data
    X = np.concatenate((X,5*np.random.random_sample((150,p))),1)
    y = data.target
    test_skdatasets_classif_MIFilter_vs_f_classif(X, y,p+k, k)


def test_make_classification_MIFilter_vs_f_classif():
    #This initially creates clusters of points normally distributed (std=1) about vertices of a 2 * class_sep-sided
    #hypercube, and assigns an equal number of clusters to each class. It introduces interdependence between
    #these features and adds various types of further noise to the data.
    print("make_classification Dataset\n")
    #X, y = datasets.make_classification(n_samples=200, n_features=20, n_informative=5, shuffle=False, random_state=0)
    #X, y = datasets.make_classification(n_samples=200, n_features=20, n_informative=3, n_redundant=2,n_repeated=0,
    #                        n_classes=8, n_clusters_per_class=1, flip_y=0.0,class_sep=10, shuffle=False, random_state=0)
    #p=20; k = 5
    X, y = datasets.make_classification(n_samples=200, n_features=3, n_informative=1, n_redundant=0,n_repeated=0,
                            n_classes=2, n_clusters_per_class=1, flip_y=0.0,class_sep=10, shuffle=False, random_state=0)
    p=3; k = 1
    test_skdatasets_classif_MIFilter_vs_f_classif(X, y,p+k, k)



if __name__ == '__main__':
    np.random.seed(0)
    test_equal_dependency_bivariate_data()
    test_correlation_examples(500)
    test_iris_classif_MIFilter_vs_f_classif()
    test_make_classification_MIFilter_vs_f_classif()
