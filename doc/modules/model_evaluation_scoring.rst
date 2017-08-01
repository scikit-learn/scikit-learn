.. currentmodule:: sklearn

.. _scoring_parameter:

The ``scoring`` parameter: defining model evaluation rules
==========================================================


Model selection and evaluation using tools, such as
:class:`model_selection.GridSearchCV` and
:func:`model_selection.cross_val_score`, take a ``scoring`` parameter that
controls what metric they apply to the estimators evaluated.

Common cases: predefined values
-------------------------------

For the most common use cases, you can designate a scorer object with the
``scoring`` parameter; the table below shows all possible values.
All scorer objects follow the convention that **higher return values are better
than lower return values**.  Thus metrics which measure the distance between
the model and the data, like :func:`metrics.mean_squared_error`, are
available as neg_mean_squared_error which return the negated value
of the metric.

==============================    =============================================     ==================================
Scoring                           Function                                          Comment
==============================    =============================================     ==================================
**Classification**
'accuracy'                        :func:`metrics.accuracy_score`
'average_precision'               :func:`metrics.average_precision_score`
'f1'                              :func:`metrics.f1_score`                          for binary targets
'f1_micro'                        :func:`metrics.f1_score`                          micro-averaged
'f1_macro'                        :func:`metrics.f1_score`                          macro-averaged
'f1_weighted'                     :func:`metrics.f1_score`                          weighted average
'f1_samples'                      :func:`metrics.f1_score`                          by multilabel sample
'neg_log_loss'                    :func:`metrics.log_loss`                          requires ``predict_proba`` support
'precision' etc.                  :func:`metrics.precision_score`                   suffixes apply as with 'f1'
'recall' etc.                     :func:`metrics.recall_score`                      suffixes apply as with 'f1'
'roc_auc'                         :func:`metrics.roc_auc_score`

**Clustering**
'adjusted_mutual_info_score'      :func:`metrics.adjusted_mutual_info_score`
'adjusted_rand_score'             :func:`metrics.adjusted_rand_score`
'completeness_score'              :func:`metrics.completeness_score`
'fowlkes_mallows_score'           :func:`metrics.fowlkes_mallows_score`
'homogeneity_score'               :func:`metrics.homogeneity_score`
'mutual_info_score'               :func:`metrics.mutual_info_score`
'normalized_mutual_info_score'    :func:`metrics.normalized_mutual_info_score`
'v_measure_score'                 :func:`metrics.v_measure_score`

**Regression**
'neg_mean_absolute_error'         :func:`metrics.mean_absolute_error`
'neg_mean_squared_error'          :func:`metrics.mean_squared_error`
'neg_mean_squared_log_error'      :func:`metrics.mean_squared_log_error`
'neg_median_absolute_error'       :func:`metrics.median_absolute_error`
'r2'                              :func:`metrics.r2_score`
==============================    =============================================     ==================================


Usage examples:

    >>> from sklearn import svm, datasets
    >>> from sklearn.model_selection import cross_val_score
    >>> iris = datasets.load_iris()
    >>> X, y = iris.data, iris.target
    >>> clf = svm.SVC(probability=True, random_state=0)
    >>> cross_val_score(clf, X, y, scoring='neg_log_loss') # doctest: +ELLIPSIS
    array([-0.07..., -0.16..., -0.06...])
    >>> model = svm.SVC()
    >>> cross_val_score(model, X, y, scoring='wrong_choice')
    Traceback (most recent call last):
    ValueError: 'wrong_choice' is not a valid scoring value. Valid options are ['accuracy', 'adjusted_mutual_info_score', 'adjusted_rand_score', 'average_precision', 'completeness_score', 'f1', 'f1_macro', 'f1_micro', 'f1_samples', 'f1_weighted', 'fowlkes_mallows_score', 'homogeneity_score', 'mutual_info_score', 'neg_log_loss', 'neg_mean_absolute_error', 'neg_mean_squared_error', 'neg_mean_squared_log_error', 'neg_median_absolute_error', 'normalized_mutual_info_score', 'precision', 'precision_macro', 'precision_micro', 'precision_samples', 'precision_weighted', 'r2', 'recall', 'recall_macro', 'recall_micro', 'recall_samples', 'recall_weighted', 'roc_auc', 'v_measure_score']

.. note::

    The values listed by the ValueError exception correspond to the functions measuring
    prediction accuracy described in the following sections.
    The scorer objects for those functions are stored in the dictionary
    ``sklearn.metrics.SCORERS``.


.. _scoring:

Defining your scoring strategy from metric functions
----------------------------------------------------

The module :mod:`sklearn.metrics` also exposes a set of simple functions
measuring a prediction error given ground truth and prediction:

- functions ending with ``_score`` return a value to
  maximize, the higher the better.

- functions ending with ``_error`` or ``_loss`` return a
  value to minimize, the lower the better.  When converting
  into a scorer object using :func:`make_scorer`, set
  the ``greater_is_better`` parameter to False (True by default; see the
  parameter description below).

Metrics available for various machine learning tasks are detailed in sections
below.

Many metrics are not given names to be used as ``scoring`` values,
sometimes because they require additional parameters, such as
:func:`fbeta_score`. In such cases, you need to generate an appropriate
scoring object.  The simplest way to generate a callable object for scoring
is by using :func:`make_scorer`. That function converts metrics
into callables that can be used for model evaluation.

One typical use case is to wrap an existing metric function from the library
with non-default values for its parameters, such as the ``beta`` parameter for
the :func:`fbeta_score` function::

    >>> from sklearn.metrics import fbeta_score, make_scorer
    >>> ftwo_scorer = make_scorer(fbeta_score, beta=2)
    >>> from sklearn.model_selection import GridSearchCV
    >>> from sklearn.svm import LinearSVC
    >>> grid = GridSearchCV(LinearSVC(), param_grid={'C': [1, 10]}, scoring=ftwo_scorer)

The second use case is to build a completely custom scorer object
from a simple python function using :func:`make_scorer`, which can
take several parameters:

* the python function you want to use (``my_custom_loss_func``
  in the example below)

* whether the python function returns a score (``greater_is_better=True``,
  the default) or a loss (``greater_is_better=False``).  If a loss, the output
  of the python function is negated by the scorer object, conforming to
  the cross validation convention that scorers return higher values for better models.

* for classification metrics only: whether the python function you provided requires continuous decision
  certainties (``needs_threshold=True``).  The default value is
  False.

* any additional parameters, such as ``beta`` or ``labels`` in :func:`f1_score`.

Here is an example of building custom scorers, and of using the
``greater_is_better`` parameter::

    >>> import numpy as np
    >>> def my_custom_loss_func(ground_truth, predictions):
    ...     diff = np.abs(ground_truth - predictions).max()
    ...     return np.log(1 + diff)
    ...
    >>> # loss_func will negate the return value of my_custom_loss_func,
    >>> #  which will be np.log(2), 0.693, given the values for ground_truth
    >>> #  and predictions defined below.
    >>> loss  = make_scorer(my_custom_loss_func, greater_is_better=False)
    >>> score = make_scorer(my_custom_loss_func, greater_is_better=True)
    >>> ground_truth = [[1], [1]]
    >>> predictions  = [0, 1]
    >>> from sklearn.dummy import DummyClassifier
    >>> clf = DummyClassifier(strategy='most_frequent', random_state=0)
    >>> clf = clf.fit(ground_truth, predictions)
    >>> loss(clf,ground_truth, predictions) # doctest: +ELLIPSIS
    -0.69...
    >>> score(clf,ground_truth, predictions) # doctest: +ELLIPSIS
    0.69...


.. _diy_scoring:

Implementing your own scoring object
------------------------------------
You can generate even more flexible model scorers by constructing your own
scoring object from scratch, without using the :func:`make_scorer` factory.
For a callable to be a scorer, it needs to meet the protocol specified by
the following two rules:

- It can be called with parameters ``(estimator, X, y)``, where ``estimator``
  is the model that should be evaluated, ``X`` is validation data, and ``y`` is
  the ground truth target for ``X`` (in the supervised case) or ``None`` (in the
  unsupervised case).

- It returns a floating point number that quantifies the
  ``estimator`` prediction quality on ``X``, with reference to ``y``.
  Again, by convention higher numbers are better, so if your scorer
  returns loss, that value should be negated.

.. _multimetric_scoring:

Using multiple metric evaluation
--------------------------------

Scikit-learn also permits evaluation of multiple metrics in ``GridSearchCV``,
``RandomizedSearchCV`` and ``cross_validate``.

There are two ways to specify multiple scoring metrics for the ``scoring``
parameter:

- As an iterable of string metrics::
      >>> scoring = ['accuracy', 'precision']

- As a ``dict`` mapping the scorer name to the scoring function::
      >>> from sklearn.metrics import accuracy_score
      >>> from sklearn.metrics import make_scorer
      >>> scoring = {'accuracy': make_scorer(accuracy_score),
      ...            'prec': 'precision'}

Note that the dict values can either be scorer functions or one of the
predefined metric strings.

Currently only those scorer functions that return a single score can be passed
inside the dict. Scorer functions that return multiple values are not
permitted and will require a wrapper to return a single metric::

    >>> from sklearn.model_selection import cross_validate
    >>> from sklearn.metrics import confusion_matrix
    >>> # A sample toy binary classification dataset
    >>> X, y = datasets.make_classification(n_classes=2, random_state=0)
    >>> svm = LinearSVC(random_state=0)
    >>> def tp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
    >>> def tn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 0]
    >>> def fp(y_true, y_pred): return confusion_matrix(y_true, y_pred)[1, 0]
    >>> def fn(y_true, y_pred): return confusion_matrix(y_true, y_pred)[0, 1]
    >>> scoring = {'tp' : make_scorer(tp), 'tn' : make_scorer(tn),
    ...            'fp' : make_scorer(fp), 'fn' : make_scorer(fn)}
    >>> cv_results = cross_validate(svm.fit(X, y), X, y, scoring=scoring)
    >>> # Getting the test set true positive scores
    >>> print(cv_results['test_tp'])          # doctest: +NORMALIZE_WHITESPACE
    [12 13 15]
    >>> # Getting the test set false negative scores
    >>> print(cv_results['test_fn'])          # doctest: +NORMALIZE_WHITESPACE
    [5 4 1]
