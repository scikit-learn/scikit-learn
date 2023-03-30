"""
================================
ROC Curve with Visualization API
================================

Scikit-learn defines a simple API for creating visualizations for machine
learning. The key features of this API is to allow for quick plotting and
visual adjustments without recalculation. In this example, we will demonstrate
how to use the visualization API by comparing ROC curves.
"""

# %%
# Load Data and Train a SVC
# -------------------------
# First, we load the wine dataset and convert it to a binary classification
# problem. Then, we train a support vector classifier on a training dataset.
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import RocCurveDisplay
from sklearn.datasets import load_wine
from sklearn.model_selection import train_test_split

X, y = load_wine(return_X_y=True)
y = y == 2

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)
svc = SVC(random_state=42)
svc.fit(X_train, y_train)

# %%
# Plotting the ROC Curve
# ----------------------
# Next, we plot the ROC curve with a single call to
# :func:`sklearn.metrics.RocCurveDisplay.from_estimator`. The returned
# `svc_disp` object allows us to continue using the already computed ROC curve
# for the SVC in future plots.
svc_disp = RocCurveDisplay.from_estimator(svc, X_test, y_test)
plt.show()

# %%
# Training a Random Forest and Plotting the ROC Curve
# ---------------------------------------------------
# We train a random forest classifier and create a plot comparing it to the SVC
# ROC curve. Notice how `svc_disp` uses
# :func:`~sklearn.metrics.RocCurveDisplay.plot` to plot the SVC ROC curve
# without recomputing the values of the roc curve itself. Furthermore, we
# pass `alpha=0.8` to the plot functions to adjust the alpha values of the
# curves.
rfc = RandomForestClassifier(n_estimators=10, random_state=42)
rfc.fit(X_train, y_train)
ax = plt.gca()
rfc_disp = RocCurveDisplay.from_estimator(rfc, X_test, y_test, ax=ax, alpha=0.8)
svc_disp.plot(ax=ax, alpha=0.8)
plt.show()

# %%
# Avoiding recomputing predictions
# --------------------------------
# The :meth:`sklearn.metrics.RocCurveDisplay.from_estimator` method will compute
# the predictions using the estimator and the data. However, one can already have
# the predictions computed. Calling
# :meth:`sklearn.metrics.RocCurveDisplay.from_estimator` will waste time by
# recomputing the predictions. To avoid this, one can use
# :meth:`sklearn.metrics.RocCurveDisplay.from_predictions`.
y_pred_svc = svc.decision_function(X_test)
y_pred_rfc = rfc.predict_proba(X_test)
_, ax = plt.subplots()
rfc_disp = RocCurveDisplay.from_predictions(
    y_test, y_pred_rfc[:, 1], ax=ax, alpha=0.8, name=rfc.__class__.__name__
)
svc_disp = RocCurveDisplay.from_predictions(
    y_test, y_pred_svc, ax=ax, alpha=0.8, name=svc.__class__.__name__
)
plt.show()

# %%
# Plotting the ROC Curve with Cross Validation
# --------------------------------------------
# A model is usually evaluated using cross-validation. The
# :meth:`~sklearn.metrics.RocCurveDisplay.from_cv_results` method allows for
# plotting the ROC curves of each cross-validation split and an aggregated
# average ROC curve. We need to provide the output of
# :func:`sklearn.model_selection.cross_validate` where we need to store both
# the indices and the estimators.
from sklearn.model_selection import cross_validate

cv_results = cross_validate(rfc, X, y, return_estimator=True, return_indices=True)
rfc_disp = RocCurveDisplay.from_cv_results(cv_results, X, y, kind="both")
plt.show()

# %%
