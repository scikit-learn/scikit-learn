"""
=====================================================================
Faces recognition example using eigenfaces and kernel approximation
=====================================================================

This example builds a classical face recognition pipeline on the "Labeled Faces
in the Wild" (LFW) dataset, a preprocessed excerpt of which is available here:
https://www.kaggle.com/datasets/jessicali9530/lfw-dataset

We reduce the dimensionality of the face images with PCA (the eigenfaces), then
approximate the RBF kernel with :class:`~sklearn.kernel_approximation.Nystroem`
and train a :class:`~sklearn.linear_model.LogisticRegression` on the resulting
features. The full chain is wrapped in a :class:`~sklearn.pipeline.Pipeline` so
that cross-validation does not leak information from the test set. The
hyperparameters are tuned with a successive halving search
(:class:`~sklearn.model_selection.HalvingRandomSearchCV`) that minimizes the
log loss. We finally evaluate the model both quantitatively, with a
classification report and one-vs-rest ROC and precision-recall curves, and
qualitatively, by displaying the predictions and the eigenfaces.

"""

# Authors: The scikit-learn developers
# SPDX-License-Identifier: BSD-3-Clause

# %%
# Loading the dataset
# -------------------
#
# We download the Labeled Faces in the Wild (LFW) dataset and load it as numpy
# arrays. Each sample is a flattened grayscale image; the target is the identity
# of the person pictured.
from sklearn.datasets import fetch_lfw_people

lfw_people = fetch_lfw_people(min_faces_per_person=70, resize=0.4)

# Introspect the image arrays to find the shapes (for plotting).
n_samples, height, width = lfw_people.images.shape

# For machine learning we use the data directly (relative pixel positions are
# ignored by this model).
X = lfw_people.data
n_features = X.shape[1]

# The label to predict is the id of the person.
y = lfw_people.target
target_names = lfw_people.target_names
n_classes = target_names.shape[0]

# %%
print("Total dataset size:")
print(f"n_samples: {n_samples}")
print(f"n_features: {n_features}")
print(f"n_classes: {n_classes}")

# %%
# Splitting the dataset
# ---------------------
#
# We hold out 25% of the data for testing. Preprocessing and model fitting are
# chained in a pipeline below so that scaling and feature extraction are learned
# only from the training folds during cross-validation.
from sklearn.model_selection import train_test_split

X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.25, random_state=42
)

# %%
# Building the model pipeline
# ---------------------------
#
# We chain preprocessing and classification in a
# :class:`~sklearn.pipeline.Pipeline`. PCA extracts eigenfaces as a compact
# representation; :class:`~sklearn.kernel_approximation.Nystroem` approximates
# the RBF feature map so that a linear
# :class:`~sklearn.linear_model.LogisticRegression` can model non-linear
# decision boundaries while scaling better than a kernel SVM. Because logistic
# regression outputs calibrated probabilities, we can tune the model by
# minimizing the log loss.
from sklearn.decomposition import PCA
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

n_components = 150

model = Pipeline(
    steps=[
        ("scaler", StandardScaler()),
        ("pca", PCA(n_components=n_components, svd_solver="randomized", whiten=True)),
        ("nystroem", Nystroem(random_state=42)),
        ("logreg", LogisticRegression(max_iter=5_000)),
    ]
)
model

# %%
# Tuning the pipeline with successive halving
# -------------------------------------------
#
# We tune the ``gamma`` and ``n_components`` of the Nystroem approximation and
# the ``C`` regularization of the logistic regression with a successive halving
# search (:class:`~sklearn.model_selection.HalvingRandomSearchCV`). The search
# minimizes the log loss (``neg_log_loss``) and screens many candidates on small
# training subsets before investing compute in the most promising ones. We set
# ``min_resources`` high enough so that PCA can always extract 150 eigenfaces,
# even in the first halving iteration.
from time import time

from scipy.stats import loguniform, randint

from sklearn.experimental import enable_halving_search_cv  # noqa: F401
from sklearn.model_selection import HalvingRandomSearchCV

print("Fitting the classifier to the training set")
t0 = time()
param_distributions = {
    "nystroem__gamma": loguniform(1e-4, 1e-1),
    "nystroem__n_components": randint(50, 200),
    "logreg__C": loguniform(1e-2, 1e2),
}
clf = HalvingRandomSearchCV(
    model,
    param_distributions,
    n_candidates=30,
    factor=3,
    min_resources=300,
    scoring="neg_log_loss",
    random_state=42,
)
clf = clf.fit(X_train, y_train)
print(f"done in {time() - t0:.3f}s")

# %%
print("Best estimator found by successive halving search:")
clf.best_estimator_

# %%
# Quantitative evaluation
# -----------------------
#
# We measure the model quality on the held-out test set with a classification
# report and, since the probabilities are well calibrated, one-vs-rest ROC and
# precision-recall curves. The pipeline handles preprocessing internally.
import matplotlib.pyplot as plt

from sklearn.metrics import classification_report
from sklearn.preprocessing import label_binarize

print("Predicting people's names on the test set")
t0 = time()
y_pred = clf.predict(X_test)
y_score = clf.predict_proba(X_test)
print(f"done in {time() - t0:.3f}s")

print(classification_report(y_test, y_pred, target_names=target_names))

# %%
# Because the problem is multiclass, we summarize the ranking quality of the
# predicted probabilities with one-vs-rest curves: each identity is in turn
# treated as the positive class against all the others. The ROC curve relates the
# true positive rate to the false positive rate, while the precision-recall curve
# is more informative when the positive class is rare, as is the case here where
# each identity is a small fraction of the test set.
from sklearn.metrics import PrecisionRecallDisplay, RocCurveDisplay

classes = list(range(n_classes))
y_onehot_test = label_binarize(y_test, classes=classes)

fig, (ax_roc, ax_pr) = plt.subplots(1, 2, figsize=(13, 6))

for class_id, name in enumerate(target_names):
    RocCurveDisplay.from_predictions(
        y_onehot_test[:, class_id],
        y_score[:, class_id],
        name=name,
        ax=ax_roc,
        plot_chance_level=(class_id == n_classes - 1),
    )
    PrecisionRecallDisplay.from_predictions(
        y_onehot_test[:, class_id],
        y_score[:, class_id],
        name=name,
        ax=ax_pr,
    )

ax_roc.set_title("One-vs-rest ROC curves")
ax_pr.set_title("One-vs-rest precision-recall curves")
plt.tight_layout()
plt.show()

# %%
# Qualitative evaluation
# ----------------------
#
# We visualize a gallery of test portraits with their predicted and true labels
# to inspect the model's mistakes at a glance.


def plot_gallery(images, titles, height, width, n_row=3, n_col=4):
    """Plot a gallery of portraits."""
    fig, axs = plt.subplots(n_row, n_col, figsize=(1.8 * n_col, 2.4 * n_row))
    fig.subplots_adjust(bottom=0, left=0.01, right=0.99, top=0.90, hspace=0.35)
    for ax, image, title in zip(axs.ravel(), images, titles):
        ax.imshow(image.reshape((height, width)), cmap=plt.cm.gray)
        ax.set_title(title, size=12)
        ax.set_xticks(())
        ax.set_yticks(())
    return fig


def make_title(y_pred, y_test, target_names, i):
    pred_name = target_names[y_pred[i]].rsplit(" ", 1)[-1]
    true_name = target_names[y_test[i]].rsplit(" ", 1)[-1]
    return f"predicted: {pred_name}\ntrue:      {true_name}"


prediction_titles = [
    make_title(y_pred, y_test, target_names, i) for i in range(y_pred.shape[0])
]

plot_gallery(X_test, prediction_titles, height, width)

# %%
# Eigenfaces gallery
# ------------------
#
# We display the most significant eigenfaces, i.e. the principal components that
# form the basis of the face representation learned by the fitted pipeline.

pca = clf.best_estimator_.named_steps["pca"]
eigenfaces = pca.components_.reshape((pca.n_components_, height, width))

eigenface_titles = [f"eigenface {i}" for i in range(eigenfaces.shape[0])]
plot_gallery(eigenfaces, eigenface_titles, height, width)

plt.show()

# %%
# Conclusion
# ----------
#
# This example walks through a classical face recognition pipeline in scikit-learn:
#
# - **Eigenfaces (PCA)** reduce the high-dimensional pixel space to a compact set
#   of uncorrelated features that capture the main variations across faces.
# - **Nystroem + LogisticRegression** approximate a non-linear RBF kernel with a
#   linear model that scales better than a kernel SVM and is tuned to minimize
#   the log loss.
# - **Pipeline** chains preprocessing and classification so that cross-validation
#   does not leak information from the test set.
# - **Quantitative and qualitative evaluation** on a held-out test set confirm
#   whether the pipeline generalizes. The one-vs-rest ROC and precision-recall
#   curves show how well the predicted probabilities rank each identity against
#   the others, independently of any single decision threshold.
#
# In practice, face recognition is often better addressed with convolutional
# neural networks, but this family of models is outside the scope of the
# scikit-learn library. Interested readers should instead try PyTorch or
# TensorFlow to implement such models.
