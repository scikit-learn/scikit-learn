.. currentmodule:: sklearn.feature_selection

.. _feature_selection:

=================
Feature selection
=================


The classes in the :mod:`sklearn.feature_selection` module can be used
for feature selection/dimensionality reduction on sample sets, either to
improve estimators' accuracy scores or to boost their performance on very
high-dimensional datasets.


.. _variance_threshold:

Removing features with low variance
===================================

:class:`VarianceThreshold` is a simple baseline approach to feature selection.
It removes all features whose variance doesn't meet some threshold.
By default, it removes all zero-variance features,
i.e. features that have the same value in all samples.

As an example, suppose that we have a dataset with boolean features,
and we want to remove all features that are either one or zero (on or off)
in more than 80% of the samples.
Boolean features are Bernoulli random variables,
and the variance of such variables is given by

.. math:: \mathrm{Var}[X] = p(1 - p)

so we can select using the threshold ``.8 * (1 - .8)``::

  >>> from sklearn.feature_selection import VarianceThreshold
  >>> X = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 1, 1], [0, 1, 0], [0, 1, 1]]
  >>> sel = VarianceThreshold(threshold=(.8 * (1 - .8)))
  >>> sel.fit_transform(X)
  array([[0, 1],
         [1, 0],
         [0, 0],
         [1, 1],
         [1, 0],
         [1, 1]])

As expected, ``VarianceThreshold`` has removed the first column,
which has a probability :math:`p = 5/6 > .8` of containing a zero.

.. _univariate_feature_selection:

Univariate feature selection
============================

Univariate feature selection works by selecting the best features based on
univariate statistical tests. It can be seen as a preprocessing step
to an estimator. Scikit-learn exposes feature selection routines
as objects that implement the ``transform`` method:

 * :class:`SelectKBest` removes all but the :math:`k` highest scoring features

 * :class:`SelectPercentile` removes all but a user-specified highest scoring
   percentage of features

 * using common univariate statistical tests for each feature:
   false positive rate :class:`SelectFpr`, false discovery rate
   :class:`SelectFdr`, or family wise error :class:`SelectFwe`.

 * :class:`GenericUnivariateSelect` allows to perform univariate feature
   selection with a configurable strategy. This allows to select the best
   univariate selection strategy with hyper-parameter search estimator.

For instance, we can perform a :math:`\chi^2` test to the samples
to retrieve only the two best features as follows:

  >>> from sklearn.datasets import load_iris
  >>> from sklearn.feature_selection import SelectKBest
  >>> from sklearn.feature_selection import chi2
  >>> X, y = load_iris(return_X_y=True)
  >>> X.shape
  (150, 4)
  >>> X_new = SelectKBest(chi2, k=2).fit_transform(X, y)
  >>> X_new.shape
  (150, 2)

These objects take as input a scoring function that returns univariate scores
and p-values (or only scores for :class:`SelectKBest` and
:class:`SelectPercentile`):

 * For regression: :func:`r_regression`, :func:`f_regression`, :func:`mutual_info_regression`

 * For classification: :func:`chi2`, :func:`f_classif`, :func:`mutual_info_classif`

The methods based on F-test estimate the degree of linear dependency between
two random variables. On the other hand, mutual information methods can capture
any kind of statistical dependency, but being nonparametric, they require more
samples for accurate estimation.

.. topic:: Feature selection with sparse data

   If you use sparse data (i.e. data represented as sparse matrices),
   :func:`chi2`, :func:`mutual_info_regression`, :func:`mutual_info_classif`
   will deal with the data without making it dense.

.. warning::

    Beware not to use a regression scoring function with a classification
    problem, you will get useless results.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_feature_selection_plot_feature_selection.py`

    * :ref:`sphx_glr_auto_examples_feature_selection_plot_f_test_vs_mi.py`

.. _rfe:

Recursive feature elimination
=============================

Given an external estimator that assigns weights to features (e.g., the
coefficients of a linear model), the goal of recursive feature elimination (:class:`RFE`)
is to select features by recursively considering smaller and smaller sets of
features. First, the estimator is trained on the initial set of features and
the importance of each feature is obtained either through any specific attribute
(such as ``coef_``, ``feature_importances_``) or callable. Then, the least important
features are pruned from current set of features. That procedure is recursively
repeated on the pruned set until the desired number of features to select is
eventually reached.

:class:`RFECV` performs RFE in a cross-validation loop to find the optimal
number of features.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_feature_selection_plot_rfe_digits.py`: A recursive feature elimination example
      showing the relevance of pixels in a digit classification task.

    * :ref:`sphx_glr_auto_examples_feature_selection_plot_rfe_with_cross_validation.py`: A recursive feature
      elimination example with automatic tuning of the number of features
      selected with cross-validation.

.. _select_from_model:

Feature selection using SelectFromModel
=======================================

:class:`SelectFromModel` is a meta-transformer that can be used alongside any
estimator that assigns importance to each feature through a specific attribute (such as
``coef_``, ``feature_importances_``) or via an `importance_getter` callable after fitting.
The features are considered unimportant and removed if the corresponding
importance of the feature values are below the provided
``threshold`` parameter. Apart from specifying the threshold numerically,
there are built-in heuristics for finding a threshold using a string argument.
Available heuristics are "mean", "median" and float multiples of these like
"0.1*mean". In combination with the `threshold` criteria, one can use the
`max_features` parameter to set a limit on the number of features to select.

For examples on how it is to be used refer to the sections below.

.. topic:: Examples

    * :ref:`sphx_glr_auto_examples_feature_selection_plot_select_from_model_diabetes.py`

.. _l1_feature_selection:

L1-based feature selection
--------------------------

.. currentmodule:: sklearn

:ref:`Linear models <linear_model>` penalized with the L1 norm have
sparse solutions: many of their estimated coefficients are zero. When the goal
is to reduce the dimensionality of the data to use with another classifier,
they can be used along with :class:`~feature_selection.SelectFromModel`
to select the non-zero coefficients. In particular, sparse estimators useful
for this purpose are the :class:`~linear_model.Lasso` for regression, and
of :class:`~linear_model.LogisticRegression` and :class:`~svm.LinearSVC`
for classification::

  >>> from sklearn.svm import LinearSVC
  >>> from sklearn.datasets import load_iris
  >>> from sklearn.feature_selection import SelectFromModel
  >>> X, y = load_iris(return_X_y=True)
  >>> X.shape
  (150, 4)
  >>> lsvc = LinearSVC(C=0.01, penalty="l1", dual=False).fit(X, y)
  >>> model = SelectFromModel(lsvc, prefit=True)
  >>> X_new = model.transform(X)
  >>> X_new.shape
  (150, 3)

With SVMs and logistic-regression, the parameter C controls the sparsity:
the smaller C the fewer features selected. With Lasso, the higher the
alpha parameter, the fewer features selected.

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_dense_vs_sparse_data.py`.

.. _compressive_sensing:

.. topic:: **L1-recovery and compressive sensing**

   For a good choice of alpha, the :ref:`lasso` can fully recover the
   exact set of non-zero variables using only few observations, provided
   certain specific conditions are met. In particular, the number of
   samples should be "sufficiently large", or L1 models will perform at
   random, where "sufficiently large" depends on the number of non-zero
   coefficients, the logarithm of the number of features, the amount of
   noise, the smallest absolute value of non-zero coefficients, and the
   structure of the design matrix X. In addition, the design matrix must
   display certain specific properties, such as not being too correlated.

   There is no general rule to select an alpha parameter for recovery of
   non-zero coefficients. It can by set by cross-validation
   (:class:`LassoCV` or :class:`LassoLarsCV`), though this may lead to
   under-penalized models: including a small number of non-relevant
   variables is not detrimental to prediction score. BIC
   (:class:`LassoLarsIC`) tends, on the opposite, to set high values of
   alpha.

   **Reference** Richard G. Baraniuk "Compressive Sensing", IEEE Signal
   Processing Magazine [120] July 2007
   http://users.isr.ist.utl.pt/~aguiar/CS_notes.pdf


Tree-based feature selection
----------------------------

Tree-based estimators (see the :mod:`sklearn.tree` module and forest
of trees in the :mod:`sklearn.ensemble` module) can be used to compute
impurity-based feature importances, which in turn can be used to discard irrelevant
features (when coupled with the :class:`~feature_selection.SelectFromModel`
meta-transformer)::

  >>> from sklearn.ensemble import ExtraTreesClassifier
  >>> from sklearn.datasets import load_iris
  >>> from sklearn.feature_selection import SelectFromModel
  >>> X, y = load_iris(return_X_y=True)
  >>> X.shape
  (150, 4)
  >>> clf = ExtraTreesClassifier(n_estimators=50)
  >>> clf = clf.fit(X, y)
  >>> clf.feature_importances_  # doctest: +SKIP
  array([ 0.04...,  0.05...,  0.4...,  0.4...])
  >>> model = SelectFromModel(clf, prefit=True)
  >>> X_new = model.transform(X)
  >>> X_new.shape               # doctest: +SKIP
  (150, 2)

.. topic:: Examples:

    * :ref:`sphx_glr_auto_examples_ensemble_plot_forest_importances.py`: example on
      synthetic data showing the recovery of the actually meaningful
      features.

    * :ref:`sphx_glr_auto_examples_ensemble_plot_forest_importances_faces.py`: example
      on face recognition data.

.. _sequential_feature_selection:

Sequential Feature Selection
============================

Sequential Feature Selection [sfs]_ (SFS) is available in the
:class:`~sklearn.feature_selection.SequentialFeatureSelector` transformer.
SFS can be either forward or backward:

Forward-SFS is a greedy procedure that iteratively finds the best new feature
to add to the set of selected features. Concretely, we initially start with
zero features and find the one feature that maximizes a cross-validated score
when an estimator is trained on this single feature. Once that first feature
is selected, we repeat the procedure by adding a new feature to the set of
selected features. The procedure stops when the desired number of selected
features is reached, as determined by the `n_features_to_select` parameter.

Backward-SFS follows the same idea but works in the opposite direction:
instead of starting with no features and greedily adding features, we start
with *all* the features and greedily *remove* features from the set. The
`direction` parameter controls whether forward or backward SFS is used.

In general, forward and backward selection do not yield equivalent results.
Also, one may be much faster than the other depending on the requested number
of selected features: if we have 10 features and ask for 7 selected features,
forward selection would need to perform 7 iterations while backward selection
would only need to perform 3.

SFS differs from :class:`~sklearn.feature_selection.RFE` and
:class:`~sklearn.feature_selection.SelectFromModel` in that it does not
require the underlying model to expose a `coef_` or `feature_importances_`
attribute. It may however be slower considering that more models need to be
evaluated, compared to the other approaches. For example in backward
selection, the iteration going from `m` features to `m - 1` features using k-fold
cross-validation requires fitting `m * k` models, while
:class:`~sklearn.feature_selection.RFE` would require only a single fit, and
:class:`~sklearn.feature_selection.SelectFromModel` always just does a single
fit and requires no iterations.

.. topic:: Examples

    * :ref:`sphx_glr_auto_examples_feature_selection_plot_select_from_model_diabetes.py`

.. topic:: References:

   .. [sfs] Ferri et al, `Comparative study of techniques for
      large-scale feature selection
      <https://citeseerx.ist.psu.edu/doc_view/pid/5fedabbb3957bbb442802e012d829ee0629a01b6>`_.

Feature selection as part of a pipeline
=======================================

Feature selection is usually used as a pre-processing step before doing
the actual learning. The recommended way to do this in scikit-learn is
to use a :class:`~pipeline.Pipeline`::

  clf = Pipeline([
    ('feature_selection', SelectFromModel(LinearSVC(penalty="l1"))),
    ('classification', RandomForestClassifier())
  ])
  clf.fit(X, y)

In this snippet we make use of a :class:`~svm.LinearSVC`
coupled with :class:`~feature_selection.SelectFromModel`
to evaluate feature importances and select the most relevant features.
Then, a :class:`~ensemble.RandomForestClassifier` is trained on the
transformed output, i.e. using only relevant features. You can perform
similar operations with the other feature selection methods and also
classifiers that provide a way to evaluate feature importances of course.
See the :class:`~pipeline.Pipeline` examples for more details.
