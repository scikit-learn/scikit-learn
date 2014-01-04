"""
========================
Plotting Learning Curves
========================

A learning curve shows the validation and training score of a learning
algorithm for varying numbers of training samples. A learning curve
shows how much we benefit from adding more training data.

In this example, the learning curve of a PassiveAggressiveClassifier
is shown for the digits dataset. Note that the training score and the
cross-validation score are both not very good. However, the of the curve
can be found in more complex datasets very often: the training score
is very high at the beginning and decreases and the cross-validation
score is very low at the beginning and increases.
"""

# Author: Alexander Fabisch <afabisch@informatik.uni-bremen.de>
#
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from sklearn.naive_bayes import GaussianNB
from sklearn.datasets import load_digits
from sklearn.learning_curve import learning_curve


estimator = GaussianNB()
digits = load_digits()
X, y = digits.data, digits.target

plt.title("Learning Curves (Passive-Aggressive Classifier on Digits)")
plt.xlabel("Training examples")
plt.ylabel("Score")

n_samples_range, train_scores, test_scores = learning_curve(
    estimator, X, y, cv=10, n_jobs=1, verbose=1)
plt.plot(n_samples_range, train_scores, label="Training score")
plt.plot(n_samples_range, test_scores, label="Cross-validation score")

plt.legend(loc="best")
plt.show()
