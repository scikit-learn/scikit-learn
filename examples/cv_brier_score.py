"""
===========================================
Cross Validation, Prediction vs Probability
===========================================

Example of cross_val_predict function to evaluate classifier output quality
using cross-validation.

This example shows the difference for the brier loss function when using
predicted labels vs when using probability estimates.

"""
print(__doc__)

from sklearn.linear_model import LogisticRegression
from sklearn.cross_validation import cross_val_predict
from sklearn.metrics import brier_score_loss
from sklearn.datasets import make_classification

X, y = make_classification(n_samples=400, n_features=40, n_classes=2)

clf = LogisticRegression()

y_pred_labels = cross_val_predict(clf, X, y, cv=10)
y_pred_probas = cross_val_predict(clf, X, y, cv=10, predict_proba=True)


brier_label = brier_score_loss(y, y_pred_labels)

# We get the first column for each item in y_pred_probas, that is
# the probability estimates for class label '1'
brier_proba = brier_score_loss(y, y_pred_probas[:,1])



print("Brier loss for predicted targets: {:.3f}".format(brier_label))
print("Brier loss for probability estimates: {:.3f}".format(brier_proba))
