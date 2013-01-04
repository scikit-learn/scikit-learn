.. _model_evaluation:

===================
Model evaluation
===================
The :mod:`sklearn.metrics` implements score functions, performance metrics
and pairwise metrics and distance computations. Those functions are usefull to
assess the performance of an estimator under a specific criterion. Note that
in many cases, the `score` method of the underlying estimator is sufficient
and appropriate.

In this module, functions named as

  * `*_score` return a scalar value to maximize: the higher the better.
  * `*_loss` return a scalar value to minimize: the lower the better

.. _regression_metrics:

Regression metrics
==================

.. currentmodule:: sklearn.metrics

Mean absolute error
-------------------
The :func:`mean_absolute_error` function allows to compute the mean absolute
error, which is a risk function corresponding to the expected value
of the absolute error loss or :math:`l1`-norm loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the mean absolute error
(MAE) estimated over :math:`n_{\text{samples}}` is given by

.. math::

  \text{MAE}(y, \hat{y}) = \frac{1}{n_{\text{samples}}} \sum_{i}^{n_{\text{samples}}} \left| y_i - \hat{y}_i \right|.

.. topic:: References:

 * `Wikipedia - Mean absolute error
   <http://en.wikipedia.org/wiki/Mean_absolute_error>`_


Mean squared error
------------------
The :func:`mean_squared_error` function allows to compute the mean square
error, which is a risk function corresponding to the expected value
of the squared error loss or quadratic loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the mean squared error
(MSE) estimated over :math:`n_{\text{samples}}` is given by

.. math::

  \text{MSE}(y, \hat{y}) = \frac{1}{n_{\text{samples}}} \sum_{i}^{n_{\text{samples}}} (y_i - \hat{y}_i)^2.

.. topic:: References:

 * `Wikipedia - Mean squared error
   <http://en.wikipedia.org/wiki/Mean_squared_error>`_


R² score, the coefficient of determination
------------------------------------------
The :func:`r2_score` function allows to compute R², the coefficient of
determination. It provides a measure of how well future samples are likely to
be predicted by the model.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the score R² estimated
over :math:`n_{\text{samples}}` is given by

.. math::

  R^2(y, \hat{y}) = 1 - \frac{\sum_{i=0}^{n_{\text{samples}}} (y_i - \hat{y}_i)^2}{\sum_{i=0}^{n_{\text{samples}}} (y_i - \bar{y})^2}

where :math:`\bar{y} =  \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}}} y_i`.

.. topic:: References:

 * `Wikipedia - Coefficient of determination
   <http://en.wikipedia.org/wiki/Coefficient_of_determination>`_


.. TODO
  Classification metrics
  ======================

  Clustering metrics
  ======================

Dummy estimators
=================

.. currentmodule:: sklearn.dummy

When doing supervised learning, a simple sanity check consists in comparing one's
estimator against simple rules of thumb.
:class:`DummyClassifier` implements three such simple strategies for classification:

- `stratified` generates randomly predictions by respecting the training
  set's class distribution,
- `most_frequent` always predicts the most frequent label in the training set,
- `uniform` generates predictions uniformly at random.

Note that with all these strategies, the `predict` method completely ignores
the input data!

To illustrate :class:`DummyClassifier`, first let's create an imbalanced
dataset::

  >>> from sklearn.datasets import load_iris
  >>> iris = load_iris()
  >>> X, y = iris.data, iris.target
  >>> y[y != 1] = -1

Next, let's compare the accuracy of `SVC` and `most_frequent`::

  >>> from sklearn.dummy import DummyClassifier
  >>> from sklearn.svm import SVC
  >>> clf = SVC(kernel='linear', C=1).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.73...
  >>> clf = DummyClassifier(strategy='most_frequent', random_state=0).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.66...

We see that `SVC` doesn't do much better than a dummy classifier. Now, let's change
the kernel::

  >>> clf = SVC(kernel='rbf', C=1).fit(X, y)
  >>> clf.score(X, y)  # doctest: +ELLIPSIS
  0.99...

We see that the accuracy was boosted to almost 100%.

More generally, when the accuracy of a classifier is too close to random classification, it
probably means that something went wrong: features are not helpful, a
hyparameter is not correctly tuned, the classifier is suffering from class
imbalance, etc...

:class:`DummyRegressor` implements a simple rule of thumb for regression:
always predict the mean of the training targets.
