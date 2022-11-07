"""
===================================================
Recursive feature elimination with cross-validation
===================================================

A Recursive Feature Elimination (RFE) example with automatic tuning of the
number of features selected with cross-validation.

"""

# %%
# Data generation
# ---------------
#
# We build a classification task using 3 informative features. The introduction
# of 2 additional redundant (i.e. correlated) features has the effect that the
# selected features vary depending on the cross-validation fold. The remaining
# features are drawn at random and discarding them should not considerably
# affect the test scores.

from sklearn.datasets import make_classification

X, y = make_classification(
    n_samples=1000,
    n_features=25,
    n_informative=3,
    n_redundant=2,
    n_repeated=0,
    n_classes=8,
    n_clusters_per_class=1,
    random_state=0,
)

# %%
# Model training and selection
# ----------------------------
#
# We create the RFE object and compute the cross-validated scores. The scoring
# strategy "accuracy" optimizes the proportion of correctly classified samples.

from sklearn.feature_selection import RFECV
from sklearn.model_selection import StratifiedKFold
from sklearn.svm import SVC

min_features_to_select = 1  # Minimum number of features to consider
svc = SVC(kernel="linear")
cv = StratifiedKFold(10)

rfecv = RFECV(
    estimator=svc,
    step=1,
    cv=cv,
    scoring="accuracy",
    min_features_to_select=min_features_to_select,
)
rfecv.fit(X, y)

print(f"Optimal number of features: {rfecv.n_features_}")

# %%
# In the present case, the model with 3 features (which corresponds to the true
# generative model) is found to be the most optimal.
#
# Plot number of features VS. cross-validation scores
# ---------------------------------------------------

import matplotlib.pyplot as plt

n_scores = len(rfecv.cv_results_["mean_test_score"])
plt.figure()
plt.xlabel("Number of features selected")
plt.ylabel("Mean test accuracy")
plt.errorbar(
    range(min_features_to_select, n_scores + min_features_to_select),
    rfecv.cv_results_["mean_test_score"],
    yerr=rfecv.cv_results_["std_test_score"],
)
plt.show()

# %%
# From the plot above one can further notice a plateau with equally good scores
# for 3 to 5 selected features, which is the result of introducing correlated
# features. Indeed, the optimal model selected by the RFE can lie within this
# range. Above 5 selected features, the variance of the cross-validated test
# errors increase and the models' statistical performances decrease.
