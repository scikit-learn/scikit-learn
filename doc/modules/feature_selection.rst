
=======================
Feature selection
=======================



Univariate feature selection
=============================

Univariate feature selection works by selecting features based on a
univariate test statistic. Although it can seen as a preprocessing step
to an estimator, `scikit.learn` exposes an object to wrap as existing
estimator with feature selection and expose a new estimator:

.. autofunction:: scikits.learn.feature_selection.univariate_selection.UnivSelection

This object takes another estimator, a scoring function that returns
univariate p values, and a selection function that selects attributes.

Feature scoring functions
--------------------------

.. warning:: 

    A common case of non-functionning for feature selection is to use a 
    regression scoring function with a classification problem.

For classification
.......................

.. autofunction:: scikits.learn.feature_selection.univariate_selection.f_classif

For regression
.................

.. autofunction:: scikits.learn.feature_selection.univariate_selection.f_regression

Feature selection functions
----------------------------

.. autofunction:: scikits.learn.feature_selection.univariate_selection.select_k_best

.. autofunction:: scikits.learn.feature_selection.univariate_selection.select_percentile

.. autofunction:: scikits.learn.feature_selection.univariate_selection.select_fpr

.. autofunction:: scikits.learn.feature_selection.univariate_selection.select_fdr

.. autofunction:: scikits.learn.feature_selection.univariate_selection.select_fwe


Examples
----------

.. literalinclude:: ../../examples/plot_feature_selection.py


