"""
=========================
One-Class SVM versus SVDD
=========================

An example comparing the One-Class SVM and SVDD models for novelty
detection.

:ref:`Support Vector Data Description (SVDD) <svm_outlier_detection>`
and :ref:`One-Class SVM <svm_outlier_detection>` are unsupervised
algorithms that learn a decision function for novelty detection, i.e
the problem of classifying new data as similar or different to the
training set.

It can be shown that the One-Class SVM and SVDD models yield identical
results in the case of a stationary kernel, like RBF, but produce different
decision functions for non-stationary kernels, e.g. polynomial. This
example demonstrates this.

Note that it is incorrect to say that the SVDD is equivalent to the
One-Class SVM: these are different models, which just happen to coincide
for a particular family of kernels.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager
from sklearn import svm

print(__doc__)

random_state = np.random.RandomState(42)

xx, yy = np.meshgrid(np.linspace(-7, 7, 501), np.linspace(-7, 7, 501))
# Generate train data
X = 0.3 * random_state.randn(100, 2)
X_train = np.r_[X + 2, X - 2]
# Generate some regular novel observations
X = 0.3 * random_state.randn(20, 2)
X_test = np.r_[X + 2, X - 2]
# Generate some abnormal novel observations
X_outliers = random_state.uniform(low=-4, high=4, size=(20, 2))

# Define the models
nu = .1
kernels = [("RBF", dict(kernel="rbf", gamma=0.1)),
           ("Poly", dict(kernel="poly", degree=2, coef0=1.0)),
           ]

for kernel_name, kernel in kernels:

    # Use low tolerance to ensure better precision of the SVM
    # optimization procedure.
    classifiers = [("OCSVM", svm.OneClassSVM(nu=nu, tol=1e-8, **kernel)),
                   ("SVDD", svm.SVDD(nu=nu, tol=1e-8, **kernel)),
                   ]

    fig = plt.figure(figsize=(12, 5))
    fig.suptitle("One-Class SVM versus SVDD "
                 "(error train, error novel regular, error novel abnormal)")

    for i, (model_name, clf) in enumerate(classifiers):
        clf.fit(X_train)

        y_pred_train = clf.predict(X_train)
        y_pred_test = clf.predict(X_test)
        y_pred_outliers = clf.predict(X_outliers)
        n_error_train = y_pred_train[y_pred_train == -1].size
        n_error_test = y_pred_test[y_pred_test == -1].size
        n_error_outliers = y_pred_outliers[y_pred_outliers == 1].size

        ax = fig.add_subplot(1, 2, i + 1)

        # plot the line, the points, and the nearest vectors to the plane
        Z = clf.decision_function(np.c_[xx.ravel(), yy.ravel()])
        Z = Z.reshape(xx.shape)

        ax.contourf(xx, yy, Z, levels=np.linspace(Z.min(), 0, 7),
                    cmap=plt.cm.PuBu, zorder=-99)
        ax.contourf(xx, yy, Z, levels=[0, Z.max()], colors='palevioletred',
                    zorder=-98)
        a = ax.contour(xx, yy, Z, levels=[0], linewidths=2, colors='darkred',
                       zorder=-97)

        s = 40
        b1 = ax.scatter(X_train[:, 0], X_train[:, 1], s=s,
                        c='white', edgecolors='k')
        b2 = ax.scatter(X_test[:, 0], X_test[:, 1], c='blueviolet', s=s)
        c = ax.scatter(X_outliers[:, 0], X_outliers[:, 1], c='gold', s=s)
        ax.axis('tight')
        ax.set_xlim((-6, 6))
        ax.set_ylim((-6, 6))

        ax.set_title("%s %s (%d/200, %d/40, %d/40)"
                     % (model_name, kernel_name, n_error_train,
                        n_error_test, n_error_outliers))

        ax.legend([a.collections[0], b1, b2, c],
                  ["learned frontier", "training observations",
                   "new regular observations", "new abnormal observations"],
                  loc="lower right",
                  prop=matplotlib.font_manager.FontProperties(size=10))

plt.show()
