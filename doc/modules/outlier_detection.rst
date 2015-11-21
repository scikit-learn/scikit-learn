.. _outlier_detection:

=============================
Novelty and Outlier Detection
=============================

.. currentmodule:: sklearn

Many applications require being able to decide whether a new observation
belongs to the same distribution as existing observations (it is an
`inlier`), or should be considered as different (it is an outlier).
Often, this ability is used to clean real data sets. Two important
distinction must be made:

:novelty detection: 
  The training data is not polluted by outliers, and we are interested in
  detecting anomalies in new observations.

:outlier detection:
  The training data contains outliers, and we need to fit the central
  mode of the training data, ignoring the deviant observations.

The scikit-learn project provides a set of machine learning tools that
can be used both for novelty or outliers detection. This strategy is
implemented with objects learning in an unsupervised way from the data::

    estimator.fit(X_train)

new observations can then be sorted as inliers or outliers with a
`predict` method::

    estimator.predict(X_test)

Inliers are labeled 1, while outliers are labeled -1.

Novelty Detection
=================

Consider a data set of :math:`n` observations from the same
distribution described by :math:`p` features.  Consider now that we
add one more observation to that data set. Is the new observation so
different from the others that we can doubt it is regular? (i.e. does
it come from the same distribution?) Or on the contrary, is it so
similar to the other that we cannot distinguish it from the original
observations? This is the question addressed by the novelty detection
tools and methods.

In general, it is about to learn a rough, close frontier delimiting
the contour of the initial observations distribution, plotted in
embedding :math:`p`-dimensional space. Then, if further observations
lay within the frontier-delimited subspace, they are considered as
coming from the same population than the initial observations. Otherwise,
if they lay outside the frontier, we can say that they are abnormal with a
given confidence in our assessment.

There are two SVM-based approaches for that purpose:

1. :class:`svm.OneClassSVM` finds a hyperplane which separates the data from
   the origin by the largest margin.
2. :class:`svm.SVDD` finds a sphere with a minimum radius which encloses
   the data.

Both methods can implicitly work in transformed high-dimensional space using
the kernel trick, the RBF kernel is used by default. :class:`svm.OneClassSVM`
provides :math:`\nu` parameter for controlling the trade off between the
margin and the number of outliers during training, namely it is an upper bound
on the fraction of outliers in a training set or probability of finding a
new, but regular, observation outside the frontier. :clss:`svm.SVDD` provides a
similar parameter :math:`C = 1 / (\nu l)`, where :math:`l` is the number of
samples, such that :math:`1/C` approximately equals the number of outliers in
a training set.

.. topic:: References:

    * Bernhard Schölkopf et al, `Estimating the support of a high-dimensional
      distribution <http://dl.acm.org/citation.cfm?id=1119749>`_, Neural
      computation 13.7 (2001): 1443-1471.
    * David M. J. Tax and Robert P. W. Duin, `Support vector data description
      <http://dl.acm.org/citation.cfm?id=960109>`_, Machine Learning,
      54(1):45-66, 2004.
      
.. topic:: Examples:

   * See :ref:`example_svm_plot_oneclass.py` for visualizing the
     frontier learned around some data by :class:`svm.OneClassSVM`.
   * See :ref:`example_svm_plot_oneclass_vs_svdd.py` to get the idea about
     the difference between the two approaches.

.. figure:: ../auto_examples/svm/images/plot_oneclass_001.png
   :target: ../auto_examples/svm/plot_oneclasse.html
   :align: center
   :scale: 75%

.. figure:: ../auto_examples/svm/images/plot_oneclass_vs_svdd_001.png
   :target: ../auto_examples/svm/plot_oneclass_vs_svdd.html
   :align: center
   :scale: 75


Outlier Detection
=================

Outlier detection is similar to novelty detection in the sense that
the goal is to separate a core of regular observations from some
polluting ones, called "outliers". Yet, in the case of outlier
detection, we don't have a clean data set representing the population
of regular observations that can be used to train any tool.


Fitting an elliptic envelope
----------------------------

One common way of performing outlier detection is to assume that the
regular data come from a known distribution (e.g. data are Gaussian
distributed). From this assumption, we generally try to define the
"shape" of the data, and can define outlying observations as
observations which stand far enough from the fit shape. 

The scikit-learn provides an object
:class:`covariance.EllipticEnvelope` that fits a robust covariance
estimate to the data, and thus fits an ellipse to the central data
points, ignoring points outside the central mode.

For instance, assuming that the inlier data are Gaussian distributed, it
will estimate the inlier location and covariance in a robust way (i.e.
whithout being influenced by outliers). The Mahalanobis distances
obtained from this estimate is used to derive a measure of outlyingness.
This strategy is illustrated below.

.. figure:: ../auto_examples/covariance/images/plot_mahalanobis_distances_001.png
   :target: ../auto_examples/covariance/plot_mahalanobis_distances.html
   :align: center
   :scale: 75%

.. topic:: Examples:

   * See :ref:`example_covariance_plot_mahalanobis_distances.py` for
     an illustration of the difference between using a standard
     (:class:`covariance.EmpiricalCovariance`) or a robust estimate
     (:class:`covariance.MinCovDet`) of location and covariance to
     assess the degree of outlyingness of an observation.

.. topic:: References:

    ..  [RD1999] Rousseeuw, P.J., Van Driessen, K. "A fast algorithm for the minimum
        covariance determinant estimator" Technometrics 41(3), 212 (1999)


Isolation Forest
----------------

One efficient way of performing outlier detection in high-dimensional datasets
is to use random forests.
The :class:`ensemble.IsolationForest` 'isolates' observations by randomly selecting
a feature and then randomly selecting a split value between the maximum and
minimum values of the selected feature.

Since recursive partitioning can be represented by a tree structure, the
number of splittings required to isolate a sample is equivalent to the path
length from the root node to the terminating node.

This path length, averaged over a forest of such random trees, is a
measure of abnormality and our decision function.

Random partitioning produces noticeably shorter paths for anomalies.
Hence, when a forest of random trees collectively produce shorter path
lengths for particular samples, they are highly likely to be anomalies.

This strategy is illustrated below.

.. figure:: ../auto_examples/ensemble/images/plot_isolation_forest_001.png
   :target: ../auto_examples/ensemble/plot_isolation_forest.html
   :align: center
   :scale: 75%

.. topic:: Examples:

   * See :ref:`example_ensemble_plot_isolation_forest.py` for
     an illustration of the use of IsolationForest.

   * See :ref:`example_covariance_plot_outlier_detection.py` for a
     comparison of :class:`ensemble.IsolationForest` with
     :class:`svm.OneClassSVM` (tuned to perform like an outlier detection
     method) and a covariance-based outlier detection with
     :class:`covariance.MinCovDet`.

.. topic:: References:

    .. [LTZ2008] Liu, Fei Tony, Ting, Kai Ming and Zhou, Zhi-Hua. "Isolation forest."
           Data Mining, 2008. ICDM'08. Eighth IEEE International Conference on.

     
Comparison of different approaches
----------------------------------

Strictly-speaking, the SVM-based methods are not designed for outlier
detection, but rather for novelty detection: its training set should not be
contaminated by outliers as it may fit them. That said, outlier detection in
high-dimension, or without any assumptions on the distribution of the inlying
data is very challenging, and a SVM-based methods give useful results in these
situations.

The examples below illustrate how the performance of the
:class:`covariance.EllipticEnvelope` degrades as the data is less and less
unimodal, and other methods become more beneficial. Note, that the parameters
of :class:`svm.OneClassSVM` and :class:`svm.SVDD` are set to achieve their
equivalence, i. e. :math:`C = 1 / (\nu l)`.

|

- For a inlier mode well-centered and elliptic all methods give approximately
  equally good results.

.. figure:: ../auto_examples/covariance/images/plot_outlier_detection_001.png
   :target: ../auto_examples/covariance/plot_outlier_detection.html
   :align: center
   :scale: 75%

- As the inlier distribution becomes bimodal,
  :class:`covariance.EllipticEnvelope` does not fit well the inliers. However,
  we can see that other methods also have difficulties to detect the two modes,
  but generally perform equally well.

.. figure:: ../auto_examples/covariance/images/plot_outlier_detection_002.png
   :target: ../auto_examples/covariance/plot_outlier_detection.html
   :align: center
   :scale: 75%

- As the inlier distribution gets strongly non-Gaussian,
  :class:`covariance.EllipticEnvelope` starts to perform inadequate. Other
  methods give a reasonable representation, with
  :class:`ensemble.IsolationForest` having the least amount of errors.

.. figure:: ../auto_examples/covariance/images/plot_outlier_detection_003.png
   :target: ../auto_examples/covariance/plot_outlier_detection.html
   :align: center
   :scale: 75%
