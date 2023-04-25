.. currentmodule:: sklearn.model_selection

.. _cutoffclassifier:

========================================================
Tuning cut-off decision threshold for classes prediction
========================================================

Classifiers are predictive models: they use statistical learning to predict
outcomes. The outcomes of a classifier takes two forms: a "soft" score for each
sample in relation to each class, and a "hard" categorical prediction (i.e.
class label). Soft predictions are obtained using :term:`predict_proba` or
:term:`decision_function` while hard predictions are obtained using
:term:`predict`.

In scikit-learn, there is a connection between soft and hard prediction. In the
case of a binary classification, hard predictions are obtained by associating
the positive class with probability value greater than 0.5 (obtained with
:term:`predict_proba`) or decision function value greater than 0 (obtained with
:term:`decision_function`). Similar rules apply a defined for other
classification problems.

While these approaches are reasonable as default behaviors, they might not be
adapted to certain use cases. For instance, in a medical context, it might be
better to predict the positive class for a lower probability than 0.5 to not
miss any patient with a disease. However, it will come at the cost of having
more false positive predictions. In some use cases, one would like to define
the "hard" score based on a "business" metric instead of a statistical metric.
