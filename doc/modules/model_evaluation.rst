.. _model_evaluation:

===================
Model evaluation
===================
The :mod:`sklearn.metrics` implements score functions, performance metrics
and pairwise metrics and distance computations. Those functions are usefull to
assess the performance of an estimator under a specific criterion. Note that
in many cases, the ``score`` method of the underlying estimator is sufficient
and appropriate.

In this module, functions named as

  * `*_score` return a scalar value to maximize: the higher the better.
  * `*_loss` return a scalar value to minimize: the lower the better

.. TODO
   Missing from ref and doc
   matthews_corrcoef
   explained_variance_score


.. _classification_metrics:

Classification metrics
======================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` implements several losses, scores and utility
function to measure classification perfomance. Some of these are restricted to
the binary classification case:

.. autosummary::
   :template: function.rst

   auc_score
   average_precision_score
   hinge_loss
   precision_recall_curve
   roc_curve


Others have been extended to the multiclass case:

.. autosummary::
   :template: function.rst

  confusion_matrix
  f1_score
  fbeta_score
  precision_recall_fscore_support
  precision_score
  zero_one_score
  zero_one

In the following sub-sections, we will describe each of those functions.

Confusion matrix
----------------
The :func:`confusion_matrix` function computes the `confusion matrix
<http://en.wikipedia.org/wiki/Confusion_matrix>`_ to evaluate
the accuracy on a classification problem.

By definition a confusion matrix :math:`cm` is such that :math:`cm[i, j]` is
equal to the number of observations known to be in group :math:`i` but
predicted to be in group :math:`j`.

Precision, recall and F-measures
--------------------------------
Several functions allow you to analyse the precision, recall and F-measures
score:

.. autosummary::
   :template: function.rst

   f1_score
   fbeta_score
   precision_recall_curve
   precision_recall_fscore_support
   precision_score

.. topic:: Example:

  * See :ref:`example_plot_precision_recall.py`
    for an example of precision-Recall metric to evaluate the quality of the
    output of a classifier with :func:`precision_recall_curve`.



Binary classification
^^^^^^^^^^^^^^^^^^^^^

In a binary classification task, the terms ''positive'' and ''negative'' refer
to the classifier's prediction and the terms ''true'' and ''false'' refer to
whether that prediction corresponds to the external judgment (sometimes known
as the ''observation''). Given these definitions, we can formulate the
following table:

+-------------------+------------------------------------------------+
|                   |    Actual class (observation)                  |
+-------------------+---------------------+--------------------------+
|   Predicted class | tp (true positive)  | fp (false positive)      |
|   (expectation)   | Correct result      | Unexpected result        |
|                   +---------------------+--------------------------+
|                   | fn (false negative) | tn (true negative)       |
|                   | Missing result      | Correct absence of result|
+-------------------+---------------------+--------------------------+

In this context, we can define the notions of precision, recall and F-measure.

The precision is intuitively the ability of the classifier not to label as
positive a sample that is negative

.. math::

   \text{precision} = \frac{tp}{tp + fp}.

The recall is intuitively the ability of the classifier to find all the
positive samples

.. math::

   \text{recall} = \frac{tp}{tp + fn}.


The :math:`F_\beta` measure can be interpreted as a weighted harmonic mean of
the precision and recall, where a :math:`F_\beta` measure reaches
its best value at 1 and worst score at 0.

.. math::

   F_\beta = (1 + \beta^2) \frac{\text{precision} \times \text{recall}}{\text{precision} + \text{recall}}

With :math:`\beta = 1`, the :math:`F_\beta` measure leads to the
:math:`F_1` measure, wheres the recall and he precsion are equally important.

.. topic:: References:

   * `Precision and recall
     <http://en.wikipedia.org/wiki/Precision_and_recall>`_
   * `F-measure <http://en.wikipedia.org/wiki/F1_score>`_

Multiclass and multilabels classification
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In multiclass and multilabels classification task, the notions of precision,
recall and F-measures can be applied to each label independantly.

Nevertheless, these notions can be further extended:

The functions :func:`f1_score`, :func:`fbeta_score`,
:func:`precision_recall_fscore_support` and `precision_score` support an
argument ``average`` which define the type of averaging that is performed:

 * ``"macro"``: average over classes (does not take imbalance into account).
 * ``"micro"``: average over instances (takes imbalance into account).
 * ``"weighted"``: average weighted by support (takes imbalance into account).
   It can result in f1 score that is not between precision and recall.
 * ``None``: no averaging is performed.

.. warning::

  Currently those functions support only the multiclass case. But the following
  definitions will be general and will treat the multilabel case.

Let's define some notations:

   * :math:`n_\text{labels}` and :math:`n_\text{samples}` denotes respectively the
     number of labels and the number of samples.
   * :math:`\texttt{precision}_j`, :math:`\texttt{recall}_j` and
     :math:`\texttt{F\_beta}_j` are respectively the precision, the recall and
     :math:`F_\beta` measure for the :math:`j`-th label;
   * :math:`tp_j`, :math:`fp_j` and :math:`fn_j` respectively the number of
     true positive, false positive, false negative for the :math:`j`-th label;
   * :math:`y_i` is the set of true label and
     :math:`\hat{y}_i` is the set of predicted for the
     :math:`i`-th sample;

The macro precision, recall and :math:`F_\beta` are averaged over all labels

.. math::

  \texttt{macro\_{}precision} = \frac{1}{n_\text{labels}} \sum_{j=0}^{n_\text{labels} - 1} \texttt{precision}_j

.. math::

  \texttt{macro\_{}recall} = \frac{1}{n_\text{labels}} \sum_{j=0}^{n_\text{labels} - 1} \texttt{recall}_j

.. math::

  \texttt{macro\_{}F\_{}beta} = \frac{1}{n_\text{labels}} \sum_{j=0}^{n_\text{labels} - 1} \texttt{F\_beta}_j


The micro precision, recall and :math:`F_\beta` are averaged over all instance

.. math::

  \texttt{micro\_{}precision} = \frac{\sum_{j=0}^{n_\text{labels} - 1} tp_j}{\sum_{j=0}^{n_\text{labels} - 1} tp_j + \sum_{j=0}^{n_\text{labels} - 1} fp_j}

.. math::

  \texttt{micro\_{}recall} = \frac{\sum_{j=0}^{n_\text{labels} - 1} tp_j}{\sum_{j=0}^{n_\text{labels} - 1} tp_j + \sum_{j=0}^{n_\text{labels} - 1} fn_j}

.. math::

  \texttt{micro\_{}F\_{}beta} = (1 + \beta^2) \frac{\texttt{micro\_{}precision} \times  \texttt{micro\_{}recall}}{\texttt{micro\_{}precision} +  \texttt{micro\_{}recall}}


The weighted precision, recall and :math:`F_\beta` are averaged weighted by
their support

.. math::

  \texttt{weighted\_{}precision}(y,\hat{y}) &= \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} \frac{|y_i \cap \hat{y}_i|}{|y_i \cup \hat{y}_i|}

.. math::

  \texttt{weighted\_{}recall}(y,\hat{y}) &= \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} \frac{|y_i \cap \hat{y}_i|}{|\hat{y}_i|}

.. math::

  \texttt{weighted\_{}F\_{}beta}(y,\hat{y}) &= \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} (1 + \beta^2)\frac{|y_i \cap \hat{y}_i|}{|y_i| + |\hat{y}_i|}



Hinge loss
----------

The :func:`hinge_loss` function computes the average
`hinge loss function <http://en.wikipedia.org/wiki/Hinge_loss>`_. The hinge
loss is used in maximal margin classification as support vector machines.

If the labels are encoded with +1 and -1,  :math:`y`: the true
value and :math:`w`, the predicted decisions as output by
``decision_function``, then the hinge loss is defined as:

.. math::

  L(y, w) = \max\left\{1 - wy, 0\right\} = \left|1 - wy\right|_+



Zero one loss
--------------
The :func:`zero_one` function computes the 0-1 classification loss over
:math:`n_{\text{samples}}`. If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the 0-1 loss :math:`L_{0-1}` is defined as:

.. math::

   L_{0-1}(y_i, \hat{y}_i) = 1(\hat{y} \not= y)

where :math:`1(x)` is the indicator function.


.. _regression_metrics:

Zero one score, the accuraccy
-----------------------------
The :func:`zero_one` function computes the
`accuracy <http://en.wikipedia.org/wiki/Accuracy_and_precision>`_, the fraction
of correct`predictions.

If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the fraction of correct predictions over :math:`n_\text{samples}` is
defined as

.. math::

   \texttt{zero\_{}one}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples}-1} 1(\hat{y} = y)

where :math:`1(x)` is the indicator function.

Regression metrics
==================

.. currentmodule:: sklearn.metrics

Mean absolute error
-------------------
The :func:`mean_absolute_error` function computes the `mean absolute
error <http://en.wikipedia.org/wiki/Mean_absolute_error>`_, which is a risk
function corresponding to the expected value of the absolute error loss or
:math:`l1`-norm loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the mean absolute error
(MAE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MAE}(y, \hat{y}) = \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}}-1} \left| y_i - \hat{y}_i \right|.


Mean squared error
------------------
The :func:`mean_squared_error` function computes the `mean square
error <http://en.wikipedia.org/wiki/Mean_squared_error>`_, which is a risk
function corresponding to the expected value of the squared error loss or
quadratic loss.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the mean squared error
(MSE) estimated over :math:`n_{\text{samples}}` is defined as

.. math::

  \text{MSE}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples} - 1} (y_i - \hat{y}_i)^2.


R² score, the coefficient of determination
------------------------------------------
The :func:`r2_score` function computes R², the `coefficient of
determination <http://en.wikipedia.org/wiki/Coefficient_of_determination>`_.
It provides a measure of how well future samples are likely to
be predicted by the model.

If :math:`\hat{y}_i` is the predicted value of the :math:`i`-th sample
and :math:`y_i` is the corresponding true value, then the score R² estimated
over :math:`n_{\text{samples}}` is defined as

.. math::

  R^2(y, \hat{y}) = 1 - \frac{\sum_{i=0}^{n_{\text{samples}} - 1} (y_i - \hat{y}_i)^2}{\sum_{i=0}^{n_\text{samples} - 1} (y_i - \bar{y})^2}

where :math:`\bar{y} =  \frac{1}{n_{\text{samples}}} \sum_{i=0}^{n_{\text{samples}}} y_i`.

.. TODO

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
