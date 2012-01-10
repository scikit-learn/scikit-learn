
.. _feature_selection:

=================
Feature selection
=================

.. currentmodule:: sklearn.feature_selection

The classes in the :mod:`sklearn.feature_selection` module can be used
for feature selection/dimensionality reduction on sample sets, either to
improve estimators' accuracy scores or to boost their performance on very
high-dimensional datasets.

Univariate feature selection
============================

Univariate feature selection works by selecting the best features based on
univariate statistical tests. It can seen as a preprocessing step
to an estimator. Scikit-Learn exposes feature selection routines
a objects that implement the `transform` method:

 * selecting the k-best features :class:`SelectKBest` 

 * setting a percentile of features to keep :class:`SelectPercentile`

 * using common univariate statistical tests for each feature:
   false positive rate :class:`SelectFpr`, false discovery rate 
   :class:`SelectFdr`, or family wise error :class:`SelectFwe`.

These objects take as input a scoring function that returns
univariate p-values:

 * For regression: :func:`f_regression`

 * For classification: :func:`chi2` or :func:`f_classif`

.. topic:: Feature selection with sparse data

   If you use sparse data (i.e. data represented as sparse matrices),
   only :func:`chi2` will deal with the data without making it dense.

.. warning::

    Beware not to use a regression scoring function with a classification
    problem, you will get useless results.

.. topic:: Examples:

    :ref:`example_plot_feature_selection.py`


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

  >>> from sklearn.svm import LinearSVC
  >>> from sklearn.datasets import load_iris
  >>> iris = load_iris()
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


Tree-based feature selection
============================

Tree-based estimators (see the :mod:`sklearn.tree` module and forest
of trees in the :mod:`sklearn.ensemble` module) can be used to compute
feature importances, which in turn can be used to discard irrelevant
features::

  >>> from sklearn.ensemble import ExtraTreesClassifier
  >>> from sklearn.datasets import load_iris
  >>> iris = load_iris()
  >>> X, y = iris.data, iris.target
  >>> X.shape
  (150, 4)
  >>> clf = ExtraTreesClassifier(compute_importances=True, random_state=0)
  >>> X_new = clf.fit(X, y).transform(X)
  >>> X_new.shape
  (150, 2)

.. topic:: Examples:

    * :ref:`example_ensemble_plot_forest_importances.py`: example on
      synthetic data showing the recovery of the actually meaningful
      features.

    * :ref:`example_ensemble_plot_forest_importances_faces.py`: example
      on face recognition data.


