"""
=============================
Recursive feature elimination
=============================

This example demonstrates how :class:`~sklearn.feature_selection.RFE` and
:class:`~sklearn.feature_selection.RFECV` and can be used to determine the
importance of individual pixels when classifying handwritten digits.
RFE is a method that recursively removes the least significant features and retrains
the model, allowing us to rank features by their importance.

"""  # noqa: E501

# %%
# Dataset
# -------
#
# We start by loading the handwritten digits dataset. This dataset consists of 8x8
# pixel images of handwritten digits. Each pixel is treated as a feature and we
# aim to determine which pixels are most relevant for the digit classification task.
# %%
import matplotlib.pyplot as plt

from sklearn.datasets import load_digits

# Load the digits dataset
digits = load_digits()
X = digits.images.reshape((len(digits.images), -1))
y = digits.target

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
# :class:`~sklearn.feature_selection.RFECV`.

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
# Visualizing Feature Importance after RFE
# ----------------------------------------
#
# RFECV and RFE provide a ranking of the features based on their importance.
# We can visualize this ranking to gain insights into which pixels
# (or features) are deemed most significant by RFECV in the digit
# classification task.

# %%
ranking = rfecv.ranking_.reshape(digits.images[0].shape)
plt.matshow(ranking, cmap=plt.cm.Blues)
plt.colorbar()
plt.title("Ranking of pixels with RFECV")
plt.show()
