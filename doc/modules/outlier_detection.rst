.. _outlier_detection:

===================================================
Novelty and Outlier Detection
===================================================

Many applications require being able to caracterize the support of the
data distribution. Usually, this ability is used to clean real data
sets that are polluted by abnormal observations (outlier detection),
or to detect new categories in a further addition of observations to a
given data set (novelty detection). The boundaries of the support can
hence be a rough frontier, such as computed by a :class:`OneClassSVM`
object, or can be a smooth estimation of the data distribution density
for each point of the embedding space. (Setting up a threshold on the
smooth distribution density leads to the definition of a rough
frontier all the same.)

Sklearn provides a set of machine learning tools that can be used both
for novelty or outliers detection.

Novelty Detection
=================

Consider a data set of :math:`n` observations from the same
distribution described by :math:`p` features.  Consider now that we
add one more observation to that data set. Is the new observation so
different from the others that we can doubt it is regular? (i.e. does
it come from the same distribution?) Or on the contrary, is it so
similar to the other that we cannot distinguish it from the original
observations? This is the question adressed by the novelty detection
tools and methods.

In general, it is about to learn a rough, close frontier delimiting
the contour of the initial observations distribution, plotted in
embedding :math:`p`-dimensional space. Then, if further observations
lay within the frontier-delimited subspace, they are considered as
coming from the same population than the initial
observations. Otherwise, if they lay outside the frontier, we can say
that they are abnormal with a given confidence in our assessment.

The One-Class SVM has been introduced in [1] for that purpose and
implemented in the `sklearn.svm` package in the :class:`OneClassSVM`
object. It requires the choice of a kernel and a scalar parameter to
define a frontier.
The RBF kernel is usually chosen although there exist no exact formula
or algorithm to set its bandwith parameter. This is the default in the
sklearn implementation. The :math:`\nu` parameter, also known as the 
margin of the One-Class SVM, corresponds to the probability of finding
a new, but regular, observation outside the frontier.

.. topic:: Examples:

   * See :ref:`example_svm_plot_oneclass.py` for vizualizing the frontier
     learned around some data by a :class:`OneClassSVM` object.

.. figure:: ../auto_examples/svm/images/plot_oneclass_1.png
   :target: ../auto_examples/svm/plot_oneclasse.html
   :align: center
   :scale: 75%


Outlier Detection
=================

Outlier detection is similar to novelty detection in the sense that
the goal is to separate a core of regular observations from some
polutting ones, called "outliers". Yet, in the case of outlier
detection, we don't have a clean data set representing the population
of regular observations that can be used to train any tool.

One comon way of performing outlier detection is to assume that the
regular data come from a known distribution (e.g. data are Gaussian
distributed). From this assumption, we generaly try to define the
"shape" of the data, and can define outlying observations as
observations which stand far enough from the fit shape. For instance,
assuming that the data are Gaussian distributed, it is possible to
estimate the data location and covariance in a robust way
(i.e. whithout being influenced by outliers) and use the Mahalanobis
distances obtained from this estimate to derive a measure of outlyingness.

.. topic:: Examples:

   * See :ref:`example_covariance_plot_mahalanobis_distances.py` for
     an illustration of the difference between using a standard
     (:class:`EmpiricalCovariance`) or a robust estimate (:class:`MinCovDet`)
     of location and covariance to assess the degree of outlyingness of an
     observation.
   * See :ref:`example_covariance_plot_outlier_detection.py` for a comparison
     of the :class:`OneClassSVM` (tuned to perform like an outlier detection
     method) and a covariance-based outlier detection with :class:`MinCovDet`.

.. figure:: ../auto_examples/covariance/images/plot_mahalanobis_distances_1.png
   :target: ../auto_examples/covariance/plot_mahalanobis_distances.html
   :align: center
   :scale: 75%

.. figure:: ../auto_examples/covariance/images/plot_outlier_detection_1.png
   :target: ../auto_examples/covariance/plot_outlier_detection.html
   :align: center
   :scale: 50%
.. figure:: ../auto_examples/covariance/images/plot_outlier_detection_2.png
   :target: ../auto_examples/covariance/plot_outlier_detection.html
   :align: center
   :scale: 50%
.. figure:: ../auto_examples/covariance/images/plot_outlier_detection_3.png
   :target: ../auto_examples/covariance/plot_outlier_detection.html
   :align: center
   :scale: 50%

