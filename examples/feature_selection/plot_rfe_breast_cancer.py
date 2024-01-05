"""
=============================
Recursive feature elimination
=============================

A Recursive Feature Elimination (RFE) example to find the most discriminative features
in a breast cancer dataset. An optimized verson of cross-validation for RFE is also
demonstrated.

"""  # noqa: E501

# %%
# Dataset
# -------
#
# We start by loading the breast cancer dataset, which has 30 features.

# %%
from sklearn.datasets import load_breast_cancer

X, y = load_breast_cancer(return_X_y=True)

# %%
# Splitting the dataset for evaluation
# ------------------------------------
#
# To assess the benefits of feature selection with
# :class:`~sklearn.feature_selection.RFE`, we need a training set for selecting
# features and training our model, and a test set for evaluation.
# We'll allocate 80% of the data for training and 20% for testing.

# %%
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42
)

# %%
# Benchmarking SVM without Feature Selection
# ------------------------------------------
#
# Before applying :class:`~sklearn.feature_selection.RFE`, let's benchmark the
# performance of a :class:`~sklearn.svm.SVC` using all features. This will give us
# a baseline accuracy to compare against.

# %%
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC

svc = SVC(kernel="linear", C=1)
svc.fit(X_train, y_train)
y_pred = svc.predict(X_test)
accuracy_all_features = accuracy_score(y_test, y_pred)

print(f"Accuracy using all {X_train.shape[1]} features: {accuracy_all_features:.4f}")

# %%
# Feature Selection with RFE
# --------------------------
#
# Now, we'll employ :class:`~sklearn.feature_selection.RFE` to select the 10 most
# discriminative features. The goal is to determine if a reduced set of important
# features can either maintain or even improve the classifier's performance.

# %%
from sklearn.feature_selection import RFE
from sklearn.pipeline import Pipeline

num_features_to_select = 10

# Create a pipeline with feature selection followed by SVM
pipe = Pipeline(
    [
        (
            "rfe",
            RFE(
                estimator=SVC(kernel="linear", C=1),
                n_features_to_select=num_features_to_select,
            ),
        ),
        ("svc", SVC(kernel="linear", C=1)),
    ]
)

# Fit pipeline to the data
pipe.fit(X_train, y_train)

# Compute predictions on test data
y_pred_rfe = pipe.predict(X_test)

# Get accuracy of model using selected features
accuracy_selected_features = accuracy_score(y_test, y_pred_rfe)

print(
    f"Accuracy using {num_features_to_select} selected features:"
    f" {accuracy_selected_features:.4f}"
)

# %%
# Feature Selection with RFECV
# ----------------------------
#
# To search over all possible numbers of features, we'll employ
# :class:`~sklearn.feature_selection.RFECV`. Through feature selection and cross
# validation, we can confidently determine the number of features for optimal
# model performance.

# %%
from sklearn.feature_selection import RFECV

# Choose estimator
estimator = SVC(kernel="linear")

min_features_to_select = 1

# Define RFECV object
rfecv = RFECV(
    estimator,
    step=1,
    min_features_to_select=min_features_to_select,
    cv=5,
    scoring="accuracy",
    n_jobs=-1,
)

# Fit to the data
rfecv.fit(X, y)

# Extract the optimal number of features from the best estimator
optimal_num_features = rfecv.n_features_

print(f"Optimal number of features: {optimal_num_features}")

# %%
# Plot number of features VS. cross-validation scores
# ---------------------------------------------------

# %%
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
plt.title("Recursive Feature Elimination \nwith correlated features")
plt.show()

# %%
# From the plot above one can further notice a plateau of equivalent scores
# (similar mean value and overlapping errorbars) for 8 to 12 selected
# features. This is the result of introducing correlated features. The optimal
# model selected by the RFE is trained on 14 selected features, less than half
# of the original dataset. The test accuracy decreases above 14 selected
# features, that is, keeping non-informative features leads to over-fitting
# and is therefore detrimental for the statistical performance of the models.
