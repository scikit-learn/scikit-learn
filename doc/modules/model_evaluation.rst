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
  * `*_error` or `*_loss` return a scalar value to minimize: the lower the
    better

.. _classification_metrics:

Classification metrics
======================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` implements several losses, scores and utility
functions to measure classification perfomance. Some of these are restricted to
the binary classification case:

.. autosummary::
   :template: function.rst

   auc_score
   average_precision_score
   hinge_loss
   matthews_corrcoef
   precision_recall_curve
   roc_curve


Others have been extended to the multiclass case:

.. autosummary::
   :template: function.rst

  accuraccy_sscore
  classification_report
  confusion_matrix
  f1_score
  fbeta_score
  precision_recall_fscore_support
  precision_score
  zero_one_loss

In the following sub-sections, we will describe each of those functions.

Accuraccy score
---------------
The :func:`  accuraccy_sscore` function computes the
`accuracy <http://en.wikipedia.org/wiki/Accuracy_and_precision>`_, the fraction
of correct`predictions.

If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the fraction of correct predictions over :math:`n_\text{samples}` is
defined as

.. math::

   \texttt{accuraccy}(y, \hat{y}) = \frac{1}{n_\text{samples}} \sum_{i=0}^{n_\text{samples}-1} 1(\hat{y} = y)

where :math:`1(x)` is the indicator function.

.. topic:: Example:
  * See :ref:`example_plot_permutation_test_for_classification.py`
    for an example of accuraccy score usage to assess with permutations the
    significance of a classification score.

Area under the curve (AUC)
--------------------------
The :func:`auc_score` function computes the AUC which is
the area under the receiver operating characteristic (ROC) curve.

For more information see
`wipedia article on AUC
<http://en.wikipedia.org/wiki/Receiver_operating_characteristic#Area_under_curve>`_
and the :ref:`roc_metrics` section.

Average precision score
-----------------------
The :func:`average_precision_score` function computes the verage precision (AP)
from prediction scores. This score corresponds to the area under the
precision-recall curve.

For more information see
`wipedia article on average precision
<http://en.wikipedia.org/wiki/Information_retrieval#Average_precision>`_
and the :ref:`precision_recall_f_measure_metrics` section.

Confusion matrix
----------------
The :func:`confusion_matrix` function computes the `confusion matrix
<http://en.wikipedia.org/wiki/Confusion_matrix>`_ to evaluate
the accuracy on a classification problem.

By definition a confusion matrix :math:`cm` is such that :math:`cm[i, j]` is
equal to the number of observations known to be in group :math:`i` but
predicted to be in group :math:`j`.

.. topic:: Example:

  * See :ref:`example_plot_confusion_matrix.py`
    for an example of confusion matrix usage to evaluate the quality of the
    output of a classifier.

  * See :ref:`example_plot_digits_classification.py`
    for an example of confusion matrix usage in the classification of the
    hand-written digits.

  * See :ref:`example_document_classification_20newsgroups.py`
    for an example of confusion matrix usage in the classification of text
    documents.


Clasification report
--------------------
The :func:`classification_report` function build a text report showing the main
 classification metrics.

.. topic:: Example:

  * See :ref:`example_plot_digits_classification.py`
    for an example of classification report usage in the classification of the
    hand-written digits.

  * See :ref:`example_document_classification_20newsgroups.py`
    for an example of classification report usage in the classification of text
    documents.

  * See :ref:`example_grid_search_digits.py`
    for an example of classification report usage in parameter estimation using
    grid search with a nested cross-validation


.. _precision_recall_f_measure_metrics:

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

  * See :ref:`example_document_classification_20newsgroups.py`
    for an example of f1 score usage with classification of text
    documents.

  * See :ref:`example_grid_search_digits.py`
    for an example of precision and recall score usage in parameter estimation
    using grid search with a nested cross-validation

  * See :ref:`example_plot_sparse_recovery.py`
    for an example of :func:`precision_recall_curve` usage in feature selection
    for sparse linear models

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

   * `Wikipedia article on precision and recall
     <http://en.wikipedia.org/wiki/Precision_and_recall>`_
   * `Wikipedia article on F-measure <http://en.wikipedia.org/wiki/F1_score>`_

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
     true positives, false positives, false negatives for the :math:`j`-th
     label;
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


Matthews correlation coefficient
--------------------------------
The :func:`matthews_corrcoef` function computes the matthew's correlation
coefficient (MCC) for binary classes (quoting the `wikipedia article on the
matthew's correlation coefficient
<http://en.wikipedia.org/wiki/Matthews_correlation_coefficient>`_)

    The Matthews correlation coefficient is used in machine learning as a
    measure of the quality of binary (two-class) classifications. It takes
    into account true and false positives and negatives and is generally
    regarded as a balanced measure which can be used even if the classes are
    of very different sizes. The MCC is in essence a correlation coefficient
    value between -1 and +1. A coefficient of +1 represents a perfect
    prediction, 0 an average random prediction and -1 an inverse prediction.
    The statistic is also known as the phi coefficient. [source: Wikipedia]

If :math:`tp`, :math:`tn`, :math:`fp` and :math:`fn` are respectively the
number of true positives, true negatives, false positives ans false negatives,
the MCC coefficient is defined as

.. math::

  MCC = \frac{tp \times tn - fp \times fn}{\sqrt{(tp + fp)(tp + fn)(tn + fp)(tn + fn)}}

.. _roc_metrics:

Receiver operating characteristic (ROC)
---------------------------------------

The function :func:`roc_curve` computes the `receiver operating characteristic
curve, or ROC curve (quoting
wikipedia) <http://en.wikipedia.org/wiki/Receiver_operating_characteristic>`_:

  A receiver operating characteristic (ROC), or simply ROC curve, is a
  graphical plot which illustrates the performance of a binary classifier
  system as its discrimination threshold is varied. It is created by plotting
  the fraction of true positives out of the positives (TPR = true positive
  rate) vs. the fraction of false positives out of the negatives (FPR = false
  positive rate), at various threshold settings. TPR is also known as
  sensitivity, and FPR is one minus the specificity or true negative rate.

The following figure shows an example of ROC curve.

.. image:: ../auto_examples/images/plot_roc_1.png
   :target: ../auto_examples/plot_roc.html
   :scale: 75
   :align: center

.. topic:: Examples:

  * See :ref:`example_plot_roc.py`
    for an example of receiver operating characteristic (ROC) metric to
    evaluate the quality of the output of a classifier.

  * See :ref:`example_plot_roc_crossval.py`
    for an example of receiver operating characteristic (ROC) metric to
    evaluate the quality of the output of a classifier using cross-validation.

  * See :ref:`example_plot_species_distribution_modeling.py`
    for an example of receiver operating characteristic (ROC) metric usage to
    model species distribution.

Zero one loss
--------------
The :func:`zero_one_loss` function computes the 0-1 classification loss over
:math:`n_{\text{samples}}`. If :math:`\hat{y}_i` is the predicted value of
the :math:`i`-th sample and :math:`y_i` is the corresponding true value,
then the 0-1 loss :math:`L_{0-1}` is defined as:

.. math::

   L_{0-1}(y_i, \hat{y}_i) = 1(\hat{y} \not= y)

where :math:`1(x)` is the indicator function.


.. topic:: Example:

  * See :ref:`example_plot_rfe_with_cross_validation.py`
    for an example of zero one loss usage to perform recursive feature
    elimination with cross-validation.


.. _regression_metrics:

Regression metrics
==================

.. currentmodule:: sklearn.metrics

The :mod:`sklearn.metrics` implements several losses, scores and utility
functions to measure regressiion perfomance. Some of those have enhanced
to treat the multioutput case: :func:`mean_absolute_error`,
:func:`mean_absolute_error` and :func:`mean_squared_error`.


Explained variance score
------------------------
The :func:`explained_variance_score` computes the `explained variance
regression score <http://en.wikipedia.org/wiki/Explained_variation>`_.

If :math:`\hat{y}` is the estimated target output
and :math:`y` is the corresponding (correct) target output, then the explained
variance is  estimated  as follow:

.. math::

  \texttt{explained\_{}variance\_{}score} = 1 - \frac{Var\{ y - \hat{y}\}}{Var\{y\}}

The best possible score is 1.0, lower values are worse.


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

.. topic:: Examples:

  * See :ref:`example_plot_gradient_boosting_regression.py`
    for an example of mean squared error usage to
    evaluate gradient boosting regression.


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

.. topic:: Examples:

  * See :ref:`example_plot_lasso_and_elasticnet.py`
    for an example of R² score usage to
    evaluate Lasso and Elastic Net on sparse signals.

Clustering metrics
======================
The :mod:`sklearn.metrics` implements several losses, scores and utility
for more information see the :ref:`clustering_evaluation` section.


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
