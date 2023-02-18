
.. _metric_threshold_curve:

Metric threshold curve
======================

.. currentmodule:: sklearn.inspection

Metric threshold curve is a model inspection technique that can be used
for any :term:`fitted` binary classification :term:`estimator`. This is
especially useful for non-linear or opaque :term:`estimators`. The metric
threshold curve is defined to be how the threshold-dependent metric behaves
when we change the decision threshold.

Let's consider the following trained binary classification model::

  >>> import matplotlib.pyplot as plt
  >>> from sklearn.datasets import make_classification
  >>> from sklearn.ensemble import RandomForestClassifier
  >>> from sklearn.model_selection import train_test_split
  >>> from sklearn.metrics import fbeta_score
  >>> from functools import partial
  ...
  >>> X, y = make_classification(
  ...     n_samples=10_000, weights=(0.95, ), random_state=42)
  ...
  >>> X_train_clf, X_test, y_train_clf, y_test = train_test_split(
  ...     X, y, random_state=42, stratify=y)
  >>> X_train_clf, X_train_thr, y_train_clf, y_train_thr = train_test_split(
  ...     X_train_clf, y_train_clf, random_state=42, stratify=y_train_clf)
  ...
  >>> model = RandomForestClassifier(random_state=42).fit(X_train_clf, y_train_clf)
  ...
  >>> fbeta_score(y_test, model.predict(X_test), beta=2)
  0.462...

Its validation performance, measured via the threshold-dependent metric f2
score, is suboptimal because of the  default threshold of 0.5. We can futher
look into the behaviour of that metric with::

  >>> from sklearn.inspection import metric_threshold_curve
  >>> predict_proba_thr = model.predict_proba(X_train_thr)[:, 1]
  ...
  >>> f2_values, thresholds = metric_threshold_curve(
  ...     y_train_thr, predict_proba_thr, partial(fbeta_score, beta=2))
  ...
  >>> best_thr = thresholds[np.argmax(f2_values)]
  >>> best_thr
  ... 0.21
  ...
  >>> new_predict_test = (model.predict_proba(X_test)[:, 1] > best_thr).astype(int)
  >>> fbeta_score(y_test, new_predict_test, beta=2)
  ... 0.719...

Note that the new choosen threshold optimizes the f2 score in the test set.
