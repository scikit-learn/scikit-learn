
.. _feature_selection:

=================
Feature selection
=================

.. currentmodule:: sklearn.feature_selection

The classes in the ``sklearn.feature_selection`` module can be used
for feature selection/dimensionality reduction on sample sets, either to
improve estimators' accuracy scores or to boost their performance on very
high-dimensional datasets.

Univariate feature selection
============================

.. currentmodule:: sklearn.feature_selection.univariate_selection

Univariate feature selection works by selecting the best features based on
univariate statistical tests. It can seen as a preprocessing step
to an estimator. The `scikit.learn` exposes feature selection routines
a objects that implement the `transform` method. The k-best features
can be selected based on:

.. autofunction:: SelectKBest

or by setting a percentile of features to keep using

.. autofunction:: SelectPercentile

or using common univariate statistical test for each feature:

.. autofunction:: SelectFpr
.. autofunction:: SelectFdr
.. autofunction:: SelectFwe

These objects take as input a scoring function that returns
univariate p-values.

.. topic:: Examples:

    :ref:`example_plot_feature_selection.py`


Feature scoring functions
-------------------------

.. warning::

    Beware not to use a regression scoring function with a classification problem.

For classification
..................

.. autofunction:: chi2
.. autofunction:: f_classif

.. topic:: Feature selection with sparse data

   If you use sparse data (i.e. data represented as sparse matrices),
   only :func:`chi2` will deal with the data without making it dense.


For regression
..............

.. autofunction:: f_regression


Recursive feature elimination
=============================

.. currentmodule:: sklearn.feature_selection.rfe

Given an external estimator that assigns weights to features (e.g., the
coefficients of a linear model), the goal of recursive feature elimination (RFE)
is to select features by recursively considering smaller and smaller sets of
features.  First, the estimator is trained on the initial set of features and
weights are assigned to each one of them. Then, features whose absolute weights
are the smallest are pruned from the current set features. That procedure is
recursively repeated on the pruned set until the desired number of features to
select is eventually reached.

.. topic:: Examples:

    * :ref:`example_plot_rfe_digits.py`: A recursive feature elimination example
      showing the relevance of pixels in a digit classification task.

    * :ref:`example_plot_rfe_with_cross_validation.py`: A recursive feature
      elimination example with automatic tuning of the number of features
      selected with cross-validation.

L1-based feature selection
==========================

.. currentmodule:: sklearn.feature_selection

Linear models penalized with the L1 norm have sparse solutions. When the goal
is to reduce the dimensionality of the data to use with another classifier, the
`transform` method of `LogisticRegression` and `LinearSVC` can be used::

  >>> from sklearn import datasets
  >>> from sklearn.svm import LinearSVC
  >>> iris = datasets.load_iris()
  >>> X, y = iris.data, iris.target
  >>> X.shape
  (150, 4)
  >>> X_new = LinearSVC(C=1, penalty="l1", dual=False).fit_transform(X, y)
  >>> X_new.shape
  (150, 2)

The parameter C controls the sparsity: the smaller the fewer features.

.. topic:: Examples:

    * :ref:`example_document_classification_20newsgroups.py`: Comparison
      of different algorithms for document classification including L1-based
      feature selection.
