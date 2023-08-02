"""
===============================
Plot classification probability
===============================

Plot the classification probability for different classifiers. We use a 3 class
dataset, and we classify it with a Support Vector classifier, L1 and L2
penalized logistic regression with either a One-Vs-Rest or multinomial setting,
and Gaussian process classification.

Linear SVC is not a probabilistic classifier by default but it has a built-in
calibration option enabled in this example (`probability=True`).

The logistic regression with One-Vs-Rest is not a multiclass classifier out of
the box. As a result it has more trouble in separating class 2 and 3 than the
other estimators.

"""

# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm

from sklearn import datasets
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score
from sklearn.svm import SVC

iris = datasets.load_iris()
X = iris.data[:, 0:2]  # we only take the first two features for visualization
y = iris.target

n_features = X.shape[1]
n_classes = len(np.unique(y))

C = 10
kernel = 1.0 * RBF([1.0, 1.0])  # for GPC

# Create different classifiers.
classifiers = {
    "L1 logistic": LogisticRegression(
        C=C, penalty="l1", solver="saga", multi_class="multinomial", max_iter=10000
    ),
    "L2 logistic (Multinomial)": LogisticRegression(
        C=C, penalty="l2", solver="saga", multi_class="multinomial", max_iter=10000
    ),
    "L2 logistic (OvR)": LogisticRegression(
        C=C, penalty="l2", solver="saga", multi_class="ovr", max_iter=10000
    ),
    "Linear SVC": SVC(kernel="linear", C=C, probability=True, random_state=0),
    "GPC": GaussianProcessClassifier(kernel),
}

n_classifiers = len(classifiers)

fig, axs = plt.subplots(nrows=n_classifiers, ncols=n_classes, figsize=(6, 14))

for index, (name, classifier) in enumerate(classifiers.items()):
    classifier.fit(X, y)

    y_pred = classifier.predict(X)
    accuracy = accuracy_score(y, y_pred)
    print("Accuracy (train) for %s: %0.1f%% " % (name, accuracy * 100))

    for k in classifier.classes_:
        disp = DecisionBoundaryDisplay.from_estimator(
            classifier,
            X,
            plot_method="pcolormesh",
            response_method="predict_proba",
            class_label=k,
            ax=axs[index, k],
        )
        axs[index, k].set(
            xticks=(), yticks=(), ylabel=name if k == 0 else None, title=f"Class #{k}"
        )
        idx = y_pred == k
        if idx.any():
            axs[index, k].scatter(
                X[idx, 0], X[idx, 1], marker="o", c="w", edgecolor="k"
            )

fig.colorbar(
    cm.ScalarMappable(norm=None, cmap="viridis"),
    ax=axs,
    orientation="horizontal",
    label="Probability",
)

plt.show()
