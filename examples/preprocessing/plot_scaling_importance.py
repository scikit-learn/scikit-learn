"""
=============================
Importance of Feature Scaling
=============================

Feature scaling through standardization (also called Z-score normalization) can
be an important preprocessing step for many machine learning algorithms. It
involves rescaling each feature such that it has a standard deviation equals to
1 and a mean equals to 0.

Many algorithms require features to be normalized, either to ease the
convergence (such as the logistic regression) or because fitting scaled data
leads to a completely different model than using non-scaled data (such as
KNeighbors models). The latter is demoed on the first part of the present
example.

On the second part of the example we show how Principle Component Analysis (PCA)
is impacted by normalization of features. To illustrate this, we compare the
principal components found using :class:`~sklearn.decomposition.PCA` on unscaled
data with those obatined when using a
:class:`~sklearn.preprocessing.StandardScaler` to scale data first.

In the last part of the example we show the effect of the normalization on the
accuracy of a model trained on PCA-reduced data.

"""

# Author: Tyler Lanigan <tylerlanigan@gmail.com>
#         Sebastian Raschka <mail@sebastianraschka.com>
#         Arturo Amor <david-arturo.amor-quiroz@inria.fr>
# License: BSD 3 clause

# %%
# Load and prepare data
# =====================
#
# The dataset used is the :ref:`wine_dataset` available at UCI. This dataset has
# continuous features that are heterogeneous in scale due to differing
# properties that they measure (e.g. alcohol content and malic acid).

from sklearn.datasets import load_wine
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

X, y = load_wine(return_X_y=True, as_frame=True)
scaler = StandardScaler()

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.30, random_state=42
)
scaled_X_train = scaler.fit_transform(X_train)

# %%
# Effect of rescaling on a KNeighbors model
# =========================================
#
# For the sake of visualizing the decision boundary of a KNeighbors models, in
# this section we select a subset of 2 features that have values with different
# orders of magnitude.
#
# Keep in mind that using a subset of the features to train the model may likely
# leave out predictive variables, resulting in a decision boundary that does not
# represent the decisions nor the statistical performance of a model trained on
# the full set of features.

import matplotlib.pyplot as plt
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.neighbors import KNeighborsClassifier


X_plot = X[["proline", "hue"]].to_numpy()
X_plot_scaled = StandardScaler().fit_transform(X_plot)
clf = KNeighborsClassifier(n_neighbors=20)


def fit_and_plot_model(X_plot, y, clf, ax):
    clf.fit(X_plot, y)
    disp = DecisionBoundaryDisplay.from_estimator(
        clf,
        X_plot,
        response_method="predict",
        alpha=0.5,
        ax=ax,
    )
    disp.ax_.scatter(X_plot[:, 0], X_plot[:, 1], c=y, s=20, edgecolor="k")
    disp.ax_.set_xlim((X_plot[:, 0].min(), X_plot[:, 0].max()))
    disp.ax_.set_ylim((X_plot[:, 1].min(), X_plot[:, 1].max()))
    return disp.ax_


fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 6))

fit_and_plot_model(X_plot, y, clf, ax1)
ax1.set_xlabel("proline")
ax1.set_ylabel("hue")
ax1.set_title("KNN w/o scaling")

fit_and_plot_model(X_plot_scaled, y, clf, ax2)
ax2.set_xlabel("scaled proline")
ax2.set_ylabel("scaled hue")
ax2.set_title("KNN with scaling")

plt.show()

# %%
# Here the desicion boundary shows that fitting scaled or non-scaled data lead
# to completely different KNeighbors models. The reason is that the variable
# "proline" has values which vary between 0 and 1,000; whereas the variable
# "hue" varies between 1 and 10. Because of this, distances between samples are
# mostly impacted by differences in values of "proline", while values of the
# "hue" will be comparatively ignored. If one uses
# :class:`~sklearn.preprocessing.StandardScaler` to normalize this database,
# both scaled values lay approximately between -3 and 3 and the neighbors
# structure will be impacted more or less equivalently by both variables.
#
# Effect of rescaling on a PCA dimensional reduction
# ==================================================
#
# Dimensional reduction using PCA consists of finding the features that maximize
# the variance. If one feature varies more than another only because of their
# respective scales, :class:`~sklearn.decomposition.PCA` would determine that
# such feature dominates the direction of the principal components.
#
# We can inspect the first principal components on the full dataset:

from sklearn.decomposition import PCA

pca = PCA(n_components=2).fit(X_train)
scaled_pca = PCA(n_components=2).fit(scaled_X_train)
X_train_transformed = pca.transform(X_train)
X_train_std_transformed = scaled_pca.transform(scaled_X_train)

print(f"\nPC 1 without scaling:\n{pca.components_[0]}")
print(f"\nPC 1 with scaling:\n{scaled_pca.components_[0]}")

# %%
# Indeed we find that the feature #13 (corresponding to the "proline" variable)
# dominates the direction of the first principal component, being about two
# orders of magnitude above the other features. This is contrasted when
# observing the first principal component for the scaled version of the data,
# where the orders of magnitude are roughly the same across all the features.
#
# We can visualize the distribution of the principal components in both cases:

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 7))

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
# Effect of rescaling on model's accuracy
# =======================================
#
# Here we show prediction accuracies in scaled and unscaled data.

from sklearn.metrics import accuracy_score
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LogisticRegression

unscaled_clf = make_pipeline(pca, LogisticRegression())
unscaled_clf.fit(X_train, y_train)
y_pred = unscaled_clf.predict(X_test)

scaled_clf = make_pipeline(scaler, pca, LogisticRegression())
scaled_clf.fit(X_train, y_train)
y_pred_scaled = scaled_clf.predict(X_test)

print("Prediction accuracy for the normal test dataset with PCA")
print(f"{accuracy_score(y_test, y_pred):.2%}")
print()
print("Prediction accuracy for the standardized test dataset with PCA")
print(f"{accuracy_score(y_test, y_pred_scaled):.2%}")

# %%
# A clear difference in prediction accuracies is observed when the data is
# scaled before PCA, as it vastly outperforms the unscaled version. This
# corresponds to the intuition obtained from the plot in the previous section,
# where the components become linearly separable when scaling before using PCA.
#
# Notice that in this case the models with scaled features perform better than
# the models with non-scaled features because all the variables are expected to
# be predictive and we rather avoid some of them being comparatively ignored.
#
# If the variables in lower scales were not predictive, one may experience a
# decrease of the performance after scaling the features: noisy features would
# contribute more to the prediction after scaling and therefore scaling would
# increase overfitting.
