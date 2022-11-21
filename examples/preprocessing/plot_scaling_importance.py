"""
=============================
Importance of Feature Scaling
=============================

Feature scaling through standardization (or Z-score normalization) can be an
important preprocessing step for many machine learning algorithms.
Standardization involves rescaling each feature such that it has a standard
deviation of one and a mean of zero.

While many algorithms (such as SVM, K-nearest neighbors, and logistic
regression) require features to be normalized, intuitively we can think of
Principle Component Analysis (PCA) as being a prime example of when
normalization is important. In PCA we are interested in the components that
maximize the variance. If one component (e.g. human height) varies less than
another (e.g. weight) because of their respective scales (meters vs. kilos), PCA
might determine that the direction of maximal variance more closely corresponds
with the 'weight' axis, if those features are not scaled. As a change in height
of one meter can be considered much more important than the change in weight of
one kilogram, this is clearly incorrect.

To illustrate this, in the present example we compare the principal components
found using :class:`~klearn.decomposition.PCA` on unscaled data with those
obatined when using a :class:`~sklearn.preprocessing.StandardScaler` to scale
data first.

"""

# Author: Tyler Lanigan <tylerlanigan@gmail.com>
#         Sebastian Raschka <mail@sebastianraschka.com>

# License: BSD 3 clause

# %%
# Load and prepare data
# =====================
#
# The dataset used is the :ref:`wine_dataset` available at UCI. This dataset has
# continuous features that are heterogeneous in scale due to differing
# properties that they measure (e.g. alcohol content and malic acid). We can
# inspect the first principal components of the dataset:

from sklearn.datasets import load_wine
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

features, target = load_wine(return_X_y=True)
scaler = StandardScaler()

X_train, X_test, y_train, y_test = train_test_split(
    features, target, test_size=0.30, random_state=42
)
scaled_X_train = scaler.fit_transform(X_train)

pca = PCA(n_components=2).fit(X_train)
scaled_pca = PCA(n_components=2).fit(scaled_X_train)

print(f"\nPC 1 without scaling:\n{pca.components_[0]}")
print(f"\nPC 1 with scaling:\n{scaled_pca.components_[0]}")

# %%
# One can notice that feature #13 dominates the direction, being two orders of
# magnitude above the other features. This is contrasted when observing
# the principal component for the scaled version of the data. In the scaled
# version, the orders of magnitude are roughly the same across all the features.
#
# Plot PCA components
# ===================

import matplotlib.pyplot as plt

X_train_transformed = pca.transform(X_train)
X_train_std_transformed = scaled_pca.transform(scaled_X_train)

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 7))

target_classes = range(0, 3)
colors = ("blue", "red", "green")
markers = ("^", "s", "o")

for target_class, color, marker in zip(target_classes, colors, markers):
    ax1.scatter(
        x=X_train_transformed[y_train == target_class, 0],
        y=X_train_transformed[y_train == target_class, 1],
        color=color,
        label=f"class {target_class}",
        alpha=0.5,
        marker=marker,
    )

    ax2.scatter(
        x=X_train_std_transformed[y_train == target_class, 0],
        y=X_train_std_transformed[y_train == target_class, 1],
        color=color,
        label=f"class {target_class}",
        alpha=0.5,
        marker=marker,
    )

ax1.set_title("Training dataset after PCA")
ax2.set_title("Standardized training dataset after PCA")

for ax in (ax1, ax2):
    ax.set_xlabel("1st principal component")
    ax.set_ylabel("2nd principal component")
    ax.legend(loc="upper right")
    ax.grid()

plt.tight_layout()

plt.show()

# %%
# Model training and scoring
# ==========================
#
# Here we show prediction accuracies in scaled and unscaled data.

from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
from sklearn.pipeline import make_pipeline

unscaled_clf = make_pipeline(pca, GaussianNB())
unscaled_clf.fit(X_train, y_train)
y_pred = unscaled_clf.predict(X_test)

scaled_clf = make_pipeline(scaler, pca, GaussianNB())
scaled_clf.fit(X_train, y_train)
y_pred_scaled = scaled_clf.predict(X_test)


print("\nPrediction accuracy for the normal test dataset with PCA")
print(f"{accuracy_score(y_test, y_pred):.2%}\n")

print("\nPrediction accuracy for the standardized test dataset with PCA")
print(f"{accuracy_score(y_test, y_pred_scaled):.2%}\n")

# %%
# A clear difference in prediction accuracies is observed when the data is
# scaled before PCA, as it vastly outperforms the unscaled version.
