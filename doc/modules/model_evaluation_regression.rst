.. currentmodule:: sklearn.metrics

.. _regression_metrics:

Regression metrics
===================


The :mod:`sklearn.metrics` module implements several loss, score, and utility
functions to measure regression performance. Some of those have been enhanced
to handle the multioutput case: :func:`mean_squared_error`,
:func:`mean_absolute_error`, :func:`explained_variance_score` and
:func:`r2_score`.


These functions have an ``multioutput`` keyword argument which specifies the
way the scores or losses for each individual target should be averaged. The
default is ``'uniform_average'``, which specifies a uniformly weighted mean
over outputs. If an ``ndarray`` of shape ``(n_outputs,)`` is passed, then its
entries are interpreted as weights and an according weighted average is
returned. If ``multioutput`` is ``'raw_values'`` is specified, then all
unaltered individual scores or losses will be returned in an array of shape
``(n_outputs,)``.


The :func:`r2_score` and :func:`explained_variance_score` accept an additional
value ``'variance_weighted'`` for the ``multioutput`` parameter. This option
leads to a weighting of each individual score by the variance of the
corresponding target variable. This setting quantifies the globally captured
unscaled variance. If the target variables are of different scale, then this
score puts more importance on well explaining the higher variance variables.
``multioutput='variance_weighted'`` is the default value for :func:`r2_score`
for backward compatibility. This will be changed to ``uniform_average`` in the
future.

.. _explained_variance_score:

Explained variance score
-------------------------

The :func:`explained_variance_score` computes the `explained variance
regression score <https://en.wikipedia.org/wiki/Explained_variation>`_.

If :math:`\hat{y}` is the estimated target output, :math:`y` the corresponding
(correct) target output, and :math:`Var` is `Variance
<https://en.wikipedia.org/wiki/Variance>`_, the square of the standard deviation,
then the explained variance is estimated as follow:

.. math::

  \texttt{explained\_{}variance}(y, \hat{y}) = 1 - \frac{Var\{ y - \hat{y}\}}{Var\{y\}}

The best possible score is 1.0, lower values are worse.

Here is a small example of usage of the :func:`explained_variance_score`
function::

    >>> from sklearn.metrics import explained_variance_score
    >>> y_true = [3, -0.5, 2, 7]
    >>> y_pred = [2.5, 0.0, 2, 8]
    >>> explained_variance_score(y_true, y_pred)  # doctest: +ELLIPSIS
    0.957...
    >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
    >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
    >>> explained_variance_score(y_true, y_pred, multioutput='raw_values')
    ... # doctest: +ELLIPSIS
    array([ 0.967...,  1.        ])
    >>> explained_variance_score(y_true, y_pred, multioutput=[0.3, 0.7])
    ... # doctest: +ELLIPSIS
    0.990...

.. _mean_absolute_error:

Mean absolute error
-------------------

The :func:`mean_absolute_error` function computes `mean absolute
error <https://en.wikipedia.org/wiki/Mean_absolute_error>`_, a risk
metric corresponding to the expected value of the absolute error loss or
:math:`l1`-norm loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample,
and :math:`y_i` is the corresponding true value, then the mean absolute error
(MAE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MAE}(y, \hat{y}) = \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}}-1} \left| y_i - \hat{y}_i \right|.

Here is a small example of usage of the :func:`mean_absolute_error` function::

  >>> from sklearn.metrics import mean_absolute_error
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> mean_absolute_error(y_true, y_pred)
  0.5
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> mean_absolute_error(y_true, y_pred)
  0.75
  >>> mean_absolute_error(y_true, y_pred, multioutput='raw_values')
  array([ 0.5,  1. ])
  >>> mean_absolute_error(y_true, y_pred, multioutput=[0.3, 0.7])
  ... # doctest: +ELLIPSIS
  0.849...

.. _mean_squared_error:

Mean squared error
-------------------

The :func:`mean_squared_error` function computes `mean square
error <https://en.wikipedia.org/wiki/Mean_squared_error>`_, a risk
metric corresponding to the expected value of the squared (quadratic) error or
loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample,
and :math:`y_i` is the corresponding true value, then the mean squared error
(MSE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MSE}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} (y_i - \hat{y}_i)^2.

Here is a small example of usage of the :func:`mean_squared_error`
function::

  >>> from sklearn.metrics import mean_squared_error
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> mean_squared_error(y_true, y_pred)
  0.375
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> mean_squared_error(y_true, y_pred)  # doctest: +ELLIPSIS
  0.7083...

.. topic:: Examples:

  * See :ref:`sphx_glr_auto_examples_ensemble_plot_gradient_boosting_regression.py`
    for an example of mean squared error usage to
    evaluate gradient boosting regression.

.. _mean_squared_log_error:

Mean squared logarithmic error
------------------------------

The :func:`mean_squared_log_error` function computes a risk metric
corresponding to the expected value of the squared logarithmic (quadratic)
error or loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample,
and :math:`y_i` is the corresponding true value, then the mean squared
logarithmic error (MSLE) estimated over :math:`n_{\text{samples}}` is
defined as

.. math::

  \text{MSLE}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} (\log_e (1 + y_i) - \log_e (1 + \hat{y}_i) )^2.

Where :math:`\log_e (x)` means the natural logarithm of :math:`x`. This metric
is best to use when targets having exponential growth, such as population
counts, average sales of a commodity over a span of years etc. Note that this
metric penalizes an under-predicted estimate greater than an over-predicted
estimate.

Here is a small example of usage of the :func:`mean_squared_log_error`
function::

  >>> from sklearn.metrics import mean_squared_log_error
  >>> y_true = [3, 5, 2.5, 7]
  >>> y_pred = [2.5, 5, 4, 8]
  >>> mean_squared_log_error(y_true, y_pred)  # doctest: +ELLIPSIS
  0.039...
  >>> y_true = [[0.5, 1], [1, 2], [7, 6]]
  >>> y_pred = [[0.5, 2], [1, 2.5], [8, 8]]
  >>> mean_squared_log_error(y_true, y_pred)  # doctest: +ELLIPSIS
  0.044...

.. _median_absolute_error:

Median absolute error
---------------------

The :func:`median_absolute_error` is particularly interesting because it is
robust to outliers. The loss is calculated by taking the median of all absolute
differences between the target and the prediction.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the median absolute error
(MedAE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MedAE}(y, \hat{y}) = \text{median}(\mid y_1 - \hat{y}_1 \mid, \ldots, \mid y_n - \hat{y}_n \mid).

The :func:`median_absolute_error` does not support multioutput.

Here is a small example of usage of the :func:`median_absolute_error`
function::

  >>> from sklearn.metrics import median_absolute_error
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> median_absolute_error(y_true, y_pred)
  0.5

.. _r2_score:

R² score, the coefficient of determination
-------------------------------------------

The :func:`r2_score` function computes R², the `coefficient of
determination <https://en.wikipedia.org/wiki/Coefficient_of_determination>`_.
It provides a measure of how well future samples are likely to
be predicted by the model. Best possible score is 1.0 and it can be negative
(because the model can be arbitrarily worse). A constant model that always
predicts the expected value of y, disregarding the input features, would get a
R^2 score of 0.0.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the score R² estimated
over :math:`n_{\text{samples}}` is defined as

.. math::

  R^2(y, \hat{y}) = 1 - \frac{\sum_{i=0}^{n_{\text{samples}} - 1} (y_i - \hat{y}_i)^2}{\sum_{i=0}^{n_\text{samples} - 1} (y_i - \bar{y})^2}

where :math:`\bar{y} =  \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}} - 1} y_i`.

Here is a small example of usage of the :func:`r2_score` function::

  >>> from sklearn.metrics import r2_score
  >>> y_true = [3, -0.5, 2, 7]
  >>> y_pred = [2.5, 0.0, 2, 8]
  >>> r2_score(y_true, y_pred)  # doctest: +ELLIPSIS
  0.948...
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> r2_score(y_true, y_pred, multioutput='variance_weighted')
  ... # doctest: +ELLIPSIS
  0.938...
  >>> y_true = [[0.5, 1], [-1, 1], [7, -6]]
  >>> y_pred = [[0, 2], [-1, 2], [8, -5]]
  >>> r2_score(y_true, y_pred, multioutput='uniform_average')
  ... # doctest: +ELLIPSIS
  0.936...
  >>> r2_score(y_true, y_pred, multioutput='raw_values')
  ... # doctest: +ELLIPSIS
  array([ 0.965...,  0.908...])
  >>> r2_score(y_true, y_pred, multioutput=[0.3, 0.7])
  ... # doctest: +ELLIPSIS
  0.925...


.. topic:: Example:

  * See :ref:`sphx_glr_auto_examples_linear_model_plot_lasso_and_elasticnet.py`
    for an example of R² score usage to
    evaluate Lasso and Elastic Net on sparse signals.
