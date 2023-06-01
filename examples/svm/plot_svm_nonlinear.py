"""
==============
Non-linear SVM
==============

Perform binary classification using non-linear SVC
with RBF, Linear, Polynomial, and Sigmoid kernels. The target to predict is the XOR of the inputs.

The color map illustrates the decision function learned by the SVC.

"""

import numpy as np
import matplotlib.pyplot as plt
from sklearn import svm

# Generate sample data
np.random.seed(0)
X = np.random.randn(300, 2)
Y = np.logical_xor(X[:, 0] > 0, X[:, 1] > 0)

# Define NuSVC models
nusvc_rbf = svm.NuSVC(kernel="rbf", gamma="auto")
nusvc_lin = svm.NuSVC(kernel="linear", gamma="auto")
nusvc_poly = svm.NuSVC(kernel="poly", gamma="auto")
nusvc_sig = svm.NuSVC(kernel="sigmoid", gamma="auto")

# Fit the models
models = [nusvc_rbf, nusvc_lin, nusvc_poly, nusvc_sig]
kernel_label = ["RBF", "Linear", "Polynomial", "Sigmoid"]

fig, axes = plt.subplots(nrows=1, ncols=len(models), figsize=(20, 5))

# Plot decision function and data for each model
for ix, model in enumerate(models):
    # Fit the model
    model.fit(X, Y)

    # Plot decision function
    xx, yy = np.meshgrid(np.linspace(-3, 3, 500), np.linspace(-3, 3, 500))
    Z = model.decision_function(np.c_[xx.ravel(), yy.ravel()])
    Z = Z.reshape(xx.shape)

    # Plot decision boundary and data points
    axes[ix].imshow(
        Z,
        interpolation="nearest",
        extent=(xx.min(), xx.max(), yy.min(), yy.max()),
        aspect="auto",
        origin="lower",
        cmap=plt.cm.PuOr_r,
    )
    contours = axes[ix].contour(
        xx, yy, Z, levels=[0], linewidths=2, linestyles="dashed"
    )
    axes[ix].scatter(
        X[:, 0], X[:, 1], s=30, c=Y, cmap=plt.cm.Paired, edgecolors="k"
    )
    axes[ix].set_xticks(())
    axes[ix].set_yticks(())
    axes[ix].set_title("{} kernel".format(kernel_label[ix]))

plt.tight_layout()
plt.show()
