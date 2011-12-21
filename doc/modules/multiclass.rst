
.. _multiclass:

====================================
Multiclass and multilabel algorithms
====================================

.. currentmodule:: sklearn.multiclass

This module implements multiclass and multilabel learning algorithms:
    - one-vs-the-rest / one-vs-all
    - one-vs-one
    - error correcting output codes

Multiclass classification means classification with more than two classes.
Multilabel classification is a different task, where a classifier is used to
predict a set of target labels for each instance; i.e., the set of target
classes is not assumed to be disjoint as in ordinary (binary or multiclass)
classification. This is also called any-of classification.

The estimators provided in this module are meta-estimators: they require a base
estimator to be provided in their constructor. For example, it is possible to
use these estimators to turn a binary classifier or a regressor into a
multiclass classifier. It is also possible to use these estimators with
multiclass estimators in the hope that their accuracy or runtime performance
improves.

.. note::

    You don't need to use these estimators unless you want to experiment with
    different multiclass strategies: all classifiers in scikit-learn support
    multiclass classification out-of-the-box. Below is a summary of the
    classifiers supported in scikit-learn grouped by the strategy used.

    - Inherently multiclass: Naive Bayes, :class:`LDA`.
    - One-Vs-One: :class:`SVC`.
    - One-Vs-All: :class:`LinearSVC`, :class:`LogisticRegression`,
      :class:`SGDClassifier`, :class:`RidgeClassifier`.

.. note::

    At the moment there are no evaluation metrics implemented for multilabel
    learnings.


One-Vs-The-Rest
===============

Also known as **one-vs-all**, this strategy consists in fitting one classifier
per class. For each classifier, the class is fitted against all the other
classes. In addition to its computational efficiency (only `n_classes`
classifiers are needed), one advantage of this approach is its
interpretability. Since each class is represented by one and one classifier
only, it is possible to gain knowledge about the class by inspecting its
corresponding classifier. This is the most commonly used strategy and is a
fair default choice. Below is an example::

  >>> from sklearn import datasets
  >>> from sklearn.multiclass import OneVsRestClassifier
  >>> from sklearn.svm import LinearSVC
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> OneVsRestClassifier(LinearSVC()).fit(X, y).predict(X)
  array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])

Multilabel learning with OvR
----------------------------

``OneVsRestClassifier`` also supports multilabel classification.
To use this feature, feed the classifier a list of tuples containing
target labels, like in the example below.


.. figure:: ../auto_examples/images/plot_multilabel_1.png
    :target: ../auto_examples/plot_multilabel.html
    :align: center
    :scale: 75%


.. topic:: Examples:

    * :ref:`example_plot_multilabel.py`


One-Vs-One
==========

This strategy consists in fitting one classifier per class pair.
At prediction time, the class which received the most votes is selected.
Since it requires to fit `n_classes * (n_classes - 1) / 2` classifiers,
this method is usually slower than one-vs-the-rest, due to its
O(n_classes^2) complexity. However, this method may be advantageous for
algorithms such as kernel algorithms which don't scale well with
`n_samples`. This is because each individual learning problem only involves
a small subset of the data whereas, with one-vs-the-rest, the complete
dataset is used `n_classes` times. Below is an example::

  >>> from sklearn import datasets
  >>> from sklearn.multiclass import OneVsOneClassifier
  >>> from sklearn.svm import LinearSVC
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> OneVsOneClassifier(LinearSVC()).fit(X, y).predict(X)
  array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])


Error-Correcting Output-Codes
=============================

Output-code based strategies are fairly different from one-vs-the-rest and
one-vs-one. With these strategies, each class is represented in a euclidean
space, where each dimension can only be 0 or 1. Another way to put it is
that each class is represented by a binary code (an array of 0 and 1). The
matrix which keeps track of the location/code of each class is called the
code book. The code size is the dimensionality of the aforementioned space.
Intuitively, each class should be represented by a code as unique as
possible and a good code book should be designed to optimize classification
accuracy. In this implementation, we simply use a randomly-generated code
book as advocated in [2] although more elaborate methods may be added in the
future.

At fitting time, one binary classifier per bit in the code book is fitted.
At prediction time, the classifiers are used to project new points in the
class space and the class closest to the points is chosen.

In this implementation, the `code_size` attribute allows the user to control
the number of classifiers which will be used. It is a percentage of the
total number of classes.

A number between 0 and 1 will require fewer classifiers than
one-vs-the-rest. In theory, log2(n_classes) / n_classes is sufficient to
represent each class unambiguously. However, in practice, it may not lead to
good accuracy since log2(n_classes) is much smaller than n_classes.

A number greater than than 1 will require more classifiers than
one-vs-the-rest. In this case, some classifiers will in theory correct for
the mistakes made by other classifiers, hence the name "error-correcting".
In practice, however, this may not happen as classifier mistakes will
typically be correlated. The error-correcting output codes have a similar
effect to bagging.

Example::

  >>> from sklearn import datasets
  >>> from sklearn.multiclass import OutputCodeClassifier
  >>> from sklearn.svm import LinearSVC
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> OutputCodeClassifier(LinearSVC(), code_size=2, random_state=0).fit(X, y).predict(X)
  array([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1,
         1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
         2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 1, 2, 2, 2,
         2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2])


.. topic:: References:

 * [1] "Solving multiclass learning problems via error-correcting ouput codes",
    Dietterich T., Bakiri G., 
    Journal of Artificial Intelligence Research 2, 
    1995.

 * [2] "The error coding method and PICTs", 
    James G., Hastie T.,
    Journal of Computational and Graphical statistics 7,
    1998.

 * [3] "The Elements of Statistical Learning",
    Hastie T., Tibshirani R., Friedman J., page 606 (second-edition)
    2008.
