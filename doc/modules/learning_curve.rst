.. _learning_curves:

=====================================================
Validation curves: plotting scores to evaluate models
=====================================================

.. currentmodule:: sklearn.learning_curve

There are several problems that can occur during learning and visualization
can help to identify these. Every algorithm has its advantages and drawbacks.
A good way of categorizing estimators are their inherent `bias and variance
<http://en.wikipedia.org/wiki/Bias-variance_dilemma>`_. Simple algorithms
(e.g. :ref:`linear_model`) usually have a high **bias**, which means they have
strong assumptions about the underlying function that will be approximated.
As a result, they need fewer training samples to learn a function that
complies with these assumptions, e.g. a linear model can approximate the
true linear function usually with only very few examples very accurately.
However, if the function to be learned does not comply with these assumptions,
the error of the model will be very high and that applies to the training error
as well as to the validation error. More training examples do not help in this
case. If both the training error and the validation error are high, we call
this **underfitting**.

An example for underfitting can be seen on the left side of the following plot.
A simple linear model can at best provide only a poor fit to samples drawn from
a sine function. With polynomial features, we can increase the complexity of
the model and decrease the bias.

.. figure:: ../auto_examples/images/plot_polynomial_regression_1.png
   :target: ../auto_examples/plot_polynomial_regression.html
   :align: center
   :scale: 50%

More complex models with a low bias can usually approximate much more complex
functions. However, too complex models usually have a high **variance**, i.e.
they are very sensitive to varying training sets. Usually we can recognize too
complex models because the training error will be very low and the validation
error will be very high. This is known as **overfitting**. In the plot you can
see an example on the right side. An overfitting estimator typically learns
the noise of the training data. To reduce this type of error, we can either
select a simpler estimator with a higher bias, regularize a complex estimator,
or if this did not help, collect more training data.

We can change the bias and variance by selecting the model, which means we
first have to choose an estimator and then we have to select an appropriate
parametrization of that estimator. In the simple one-dimensional problem that
we have seen in the example, this seems to be easy. However, in
high-dimensional spaces, models can become very difficult to visualize. For
this reason, it is often helpful to use the tools described below.

.. topic:: Examples:

   * :ref:`example_plot_polynomial_regression.py`
   * :ref:`example_plot_validation_curve.py`
   * :ref:`example_plot_learning_curve.py`


.. _validation_curve:

Validation curve
================

To validate a model we need a scoring function (see :ref:`model_evaluation`),
for example accuracy for classifiers. The proper way of choosing multiple
hyperparameters of an estimator are of course grid search or similar methods
(see :ref:`grid_search`) that select the hyperparameter with the maximum score
on a validation set or multiple validation sets. Note that if we optimized
the hyperparameters based on a validation score the validation score is biased
and not a good estimate of the generalization any longer. To get a proper
estimate of the generalization we have to compute the score on another test
set.

However, it is sometimes helpful to plot the influence of a single
hyperparameter on the training score and the validation score to find out
whether the estimator is overfitting or underfitting for some hyperparameter
values.

The function :func:`validation_curve` can help in this case::

  >>> from sklearn.learning_curve import validation_curve
  >>> from sklearn.linear_model import Ridge

  >>> train_scores, valid_scores = validation_curve(Ridge(), X, y, "alpha",
  ...                                               np.logspace(-7, 3, 3))
  >>> train_scores           # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
  array([ 0.931...,  0.931...,  0.452...])
  >>> valid_scores           # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
  array([ 0.923...,  0.923...,  0.433...])

If the training score and the validation score are both low, the estimator will
be underfitting. If the training score is high and the validation score is low,
the estimator is overfitting and otherwise it is working very well. A low
training score and a high validation score is usually not possible. All three
cases can be found in the plot below where we vary the parameter
:math:`\gamma` of an SVM on the digits dataset.

.. figure:: ../auto_examples/images/plot_validation_curve_1.png
   :target: ../auto_examples/plot_validation_curve.html
   :align: center
   :scale: 50%


.. _learning_curve:

Learning curve
==============

A learning curve shows the validation and training score of an estimator
for varying numbers of training samples. It is a tool to find out how much
we benefit from adding more training data and whether the estimator suffers
more from a variance error or a bias error. If both the validation score and
the training score converge to a value that is too low with increasing
size of the training set, we will not benefit much from more training data.
In the following plot you can see an example: naive Bayes roughly converges
to a low score.

.. figure:: ../auto_examples/images/plot_learning_curve_1.png
   :target: ../auto_examples/plot_learning_curve.html
   :align: center
   :scale: 50%

We will probably have to use an estimator or a parametrization of the
current estimator that can learn more complex concepts (i.e. has a lower
bias). If the training score is much greater than the validation score for
the maximum number of training samples, adding more training samples will
most likely increase generalization. In the following plot you can see that
the SVM could benefit from more training examples.

.. figure:: ../auto_examples/images/plot_learning_curve_2.png
   :target: ../auto_examples/plot_learning_curve.html
   :align: center
   :scale: 50%

We can use the function :func:`learning_curve` to generate the values
that are required to plot such a learning curve (number of samples
that have been used, the average scores on the training sets and the
average scores on the validation sets)::

  >>> import numpy as np
  >>> from sklearn.learning_curve import learning_curve
  >>> from sklearn.datasets import load_iris
  >>> from sklearn.svm import SVC

  >>> np.random.seed(0)
  >>> iris = load_iris()
  >>> X, y = iris.data, iris.target
  >>> indices = np.arange(y.shape[0])
  >>> np.random.shuffle(indices)
  >>> X, y = X[indices], y[indices]

  >>> train_sizes, train_scores, valid_scores = learning_curve(
  ...     SVC(kernel='linear'), X, y, train_sizes=[50, 80, 110], cv=5)
  >>> train_sizes            # doctest: +NORMALIZE_WHITESPACE
  array([ 50, 80, 110])
  >>> train_scores           # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
  array([ 0.98 , 0.99 , 0.987...])
  >>> valid_scores           # doctest: +ELLIPSIS, +NORMALIZE_WHITESPACE
  array([ 0.98 , 0.986..., 0.986...])

