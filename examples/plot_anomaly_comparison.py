# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
#         Albert Thomas <albert.thomas@telecom-paristech.fr> License: BSD 3 clause

""" ============================================================================
Comparing anomaly detection algorithms for outlier detection on toy datasets
============================================================================
This example shows characteristics of different anomaly detection algorithms on 2D datasets. Datasets contain one or two modes (regions of high density) to illustrate the ability of algorithms to cope with multimodal data.
"""
print(__doc__)

import time
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from sklearn import svm
from sklearn.datasets import make_moons, make_blobs
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor

###############################################################################
#For each dataset, 15% of samples are generated as random uniform noise. This proportion is the value 
#given to the nu parameter of the OneClassSVM and the contamination parameter of the other outlier 
#detection algorithms. Decision boundaries between inliers and outliers are displayed in black except 
#for Local Outlier Factor (LOF) as it has no predict method to be applied on new data when it is used 
#for outlier detection.
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'

# Example settings
n_samples = 300
outliers_fraction = 0.15
n_outliers = int(outliers_fraction * n_samples)
n_inliers = n_samples - n_outliers

#Define Anomaly Detection Algorithms
robust_covariance = EllipticEnvelope(contamination=outliers_fraction)
oneclass_svm = svm.OneClassSVM(nu=outliers_fraction,kernel="rbf", gamma=0.1)
isolation_forest = IsolationForest(contamination=outliers_fraction, random_state=42)
local_outlier_factor = LocalOutlierFactor(n_neighbors=35, contamination=outliers_fraction)

#Algorithm Names
robust_covariance_name = "Robust covariance"
oneclass_svm_name = "One-Class SVM"
isolation_forest_name = "Isolation Forest"
local_outlier_factor_name = "Local Outlier Factor"

# Define datasets
blobs_params = dict(random_state=0, n_samples=n_inliers, n_features=2)
datasets = [
    make_blobs(centers=[[0, 0], [0, 0]], cluster_std=0.5,
               **blobs_params)[0],
    make_blobs(centers=[[2, 2], [-2, -2]], cluster_std=[0.5, 0.5],
               **blobs_params)[0],
    make_blobs(centers=[[2, 2], [-2, -2]], cluster_std=[1.5, .3],
               **blobs_params)[0],
    4. * (make_moons(n_samples=n_samples, noise=.05, random_state=0)[0] -
          np.array([0.5, 0.25])),
    14. * (np.random.RandomState(42).rand(n_samples, 2) - 0.5)]

def output_graph(datasets, algorithm, name):
    # Compare given classifiers under given settings
    xx, yy = np.meshgrid(np.linspace(-7, 7, 150),
                         np.linspace(-7, 7, 150))
    plt.figure(figsize=(8, 12.5))
    plt.subplots_adjust(left=.02, right=.98, bottom=.001, top=.96, wspace=.05,
                        hspace=.01)
    plot_num = 1
    rng = np.random.RandomState(42)
    for i_dataset, X in enumerate(datasets):
        # Add outliers
        X = np.concatenate([X, rng.uniform(low=-6, high=6,
                           size=(n_outliers, 2))], axis=0)
        t0 = time.time()
        algorithm.fit(X)
        t1 = time.time()
        plt.subplot(1, 5, plot_num)
        if i_dataset == 0:
            plt.title(name, size=18)
        # fit the data and tag outliers
        if name == "Local Outlier Factor":
            y_pred = algorithm.fit_predict(X)
        else:
            y_pred = algorithm.fit(X).predict(X)
        # plot the levels lines and the points
        if name != "Local Outlier Factor": # LOF does not implement predict
            Z = algorithm.predict(np.c_[xx.ravel(), yy.ravel()])
            Z = Z.reshape(xx.shape)
            plt.contour(xx, yy, Z, levels=[0], linewidths=2, colors='black')
        colors = np.array(['#377eb8', '#ff7f00'])
        plt.scatter(X[:, 0], X[:, 1], s=10, color=colors[(y_pred + 1) // 2])
        plt.xlim(-7, 7)
        plt.ylim(-7, 7)
        plt.xticks(())
        plt.yticks(())
        plt.text(.99, .01, ('%.2fs' % (t1 - t0)).lstrip('0'),
                 transform=plt.gca().transAxes, size=15,
                 horizontalalignment='right')
        plot_num += 1
    figure = plt.show()
    return figure

###############################################################################
#The :class:`sklearn.svm.OneClassSVM` is known to be sensitive to outliers and thus does not perform 
#very well for outlier detection. This estimator is best suited for novelty detection when the training 
#set is not contaminated by outliers. That said, outlier detection in high-dimension, or without any 
#assumptions on the distribution of the inlying data is very challenging, and a One-class SVM might give 
#useful results in these situations depending on the value of its hyperparameters.
output_graph(datasets, oneclass_svm, oneclass_svm_name)

###############################################################################
#:class:`sklearn.covariance.EllipticEnvelope` assumes the data is Gaussian and
#learns an ellipse. It thus degrades when the data is not unimodal. Notice however that this estimator 
#is robust to outliers.
output_graph(datasets, robust_covariance, robust_covariance_name)

###############################################################################
#:class:`sklearn.ensemble.IsolationForest` and class:`sklearn.neighbors.LocalOutlierFactor` seem to 
#:perform reasonably well
#for multi-modal data sets."""
output_graph(datasets,isolation_forest, isolation_forest_name)

###############################################################################
#The advantage of
#:class:`sklearn.neighbors.LocalOutlierFactor` over the other estimators is
#shown for the third data set, where the two modes have different densities. This advantage is explained 
#by the local aspect of LOF, meaning that it only compares the score of abnormality of one sample with 
#the scores of its neighbors.
output_graph(datasets, local_outlier_factor, local_outlier_factor_name)

###############################################################################
#Finally, for the last data set, it is hard to say that one sample is more abnormal than another sample 
#as they are uniformly distributed in a hypercube. Except for the :class:`sklearn.svm.OneClassSVM` which 
#overfits a little, all estimators present decent solutions for this situation. In such a samples.
#
#While these examples give some intuition about the algorithms, this intuition might not apply to very 
#high dimensional data.
#
#Finally, note that parameters of the models have been here handpicked but that in practice they need to 
#be adjusted. In the absence of labelled data,
#the problem is completely unsupervised so model selection can be a challenge.
