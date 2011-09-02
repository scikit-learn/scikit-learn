
.. _feature_selection_doc:

=================
Feature selection
=================

The classes in the ``sklearn.feature_selection`` module can be used
for feature selection/dimensionality reduction on sample sets, either to
improve estimators' accuracy scores or to boost their performance on very
high-dimensional datasets.

.. currentmodule:: sklearn.feature_selection

Univariate feature selection
============================

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


