.. currentmodule:: sklearn.model_selection

.. _prediction_tuning:

================================================
Tuning of the decision threshold of an estimator
================================================

The real-valued decision functions, i.e. `decision_function` and
`predict_proba`, of machine-learning classifiers carry the inherited biases of
the fitted model; e.g, in a class imbalanced setting, a classifier
will naturally lean toward the most frequent class. In some other cases, the
generic objective function used to train a model is generally unaware of the
evaluation criteria used to evaluate the model; e.g., one might want to
penalized differently a false-positive and false-negative ---it will be less
detrimental to show an image without a cancer (i.e., false-positive) to a
radiologist than hidding one with a cancer (i.e, false-negtative).

In a binary classification scenario, the hard-prediction, i.e. `predict`, for a
classifier most commonly use the `predict_proba` and apply a decision threshold
at 0.5 to output a positive or negative label. Thus, this hard-prediction
suffers from the same drawbacks than the one raised in the above paragraph.

Post-tuning of the decision threshold
=====================================

:class:`CutoffClassifier` allows for post-tuning the decision threshold using
either `decision_function` or `predict_proba` and an objective metric for which
we want our threshold to be optimized for.

Fine-tune using a single objective metric
-----------------------------------------

:class:`CutoffClassifier` accepts