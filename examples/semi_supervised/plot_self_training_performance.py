"""
===================================================
Self-training: Comparing performance
===================================================
This example demonstrates the performance of the SelfTrainingClassifier.

The digits dataset is loaded, and a SVC classifier is created. Then, a
SelfTrainingClassifier is initialised, using the same SVC as its
base estimator.

The dataset contains 1797 data points, and the SelfTrainingClassifier is
trained using all 1797 data points, of which some are unlabeled. The normal SVC
is trained using only the labeled data points.

The graph shows that the SelfTrainingClassifier outperforms the normal SVC
when only few labeled data points are available.

This example extends the
:ref:`sphx_glr_auto_examples_classification_plot_digits_classification.py`
example.
"""
# Authors: Oliver Rausch    <rauscho@ethz.ch>
#          Patrice Becker   <beckerp@ethz.ch>
# License: BSD 3 clause
print(__doc__)
import numpy as np
import matplotlib.pyplot as plt
from sklearn.semi_supervised.self_training import SelfTrainingClassifier
from sklearn.utils import shuffle
from sklearn.svm import SVC
from sklearn.datasets import load_digits
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import accuracy_score

cv_folds = 8
x_values = np.array([250, 500, 750, 1000])
supervised_scores = np.empty((len(x_values), cv_folds))
self_training_scores = np.empty_like(supervised_scores)

clf = SVC(probability=True, gamma=0.001)
self_training_clf = SelfTrainingClassifier(clf, max_iter=10, threshold=0.8)

X, y = load_digits(return_X_y=True)
X, y = shuffle(X, y, random_state=42)
y_true = y.copy()

for i in range(3, -1, -1):
    lim = x_values[i]
    y[lim:] = -1

    # Perform manual cross validation
    skfolds = StratifiedKFold(n_splits=cv_folds, random_state=42)
    j = 0
    for train_index, test_index in skfolds.split(X, y):
        X_train = X[train_index]
        y_train = y[train_index]
        X_test = X[test_index]
        y_test = y[test_index]
        y_test_true = y_true[test_index]

        X_train_filtered = X_train[y_train != -1]
        y_train_filtered = y_train[y_train != -1]

        # Train the supervised SVC using only data that has labels
        clf.fit(X_train_filtered, y_train_filtered)
        y_pred = clf.predict(X_test)
        supervised_scores[i, j] = accuracy_score(y_test_true, y_pred)

        # Train the semi-supervised SVC using all available data (in the fold)
        self_training_clf.fit(X_train, y_train)
        y_pred = self_training_clf.predict(X_test)
        self_training_scores[i, j] = accuracy_score(y_test_true, y_pred)
        j = j + 1

# Plot the results
plt.figure(1)
plt.errorbar(x_values, supervised_scores.mean(axis=1),
             yerr=supervised_scores.std(axis=1),
             label='Supervised (SVC)', capsize=2)
plt.errorbar(x_values, self_training_scores.mean(axis=1),
             yerr=self_training_scores.std(axis=1),
             label='SelfTrainingClassifier using SVC', capsize=2)
plt.legend()
plt.ylabel("Accuracy")
plt.title("Comparison of classifiers on limited labeled data")
plt.xlabel("Amount of labeled samples")
plt.show()
