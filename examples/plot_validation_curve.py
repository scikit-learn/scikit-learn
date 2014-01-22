"""
==========================
Plotting Validation Curves
==========================

In this plot you can see the training scores and validation scores of an SVM
for different values of the kernel parameter gamma. For very low values of
gamma, you can see that both the training score and the validation score are
low. This is called underfitting. Medium values of gamma will result in high
values for both scores, i.e. the classifier is perfoming fairly well. If gamma
is too high, the classifier will overfit, which means that the training score
is good but the validation score is poor.
"""
print(__doc__)

import matplotlib.pyplot as plt
import numpy as np
from sklearn.datasets import load_digits
from sklearn.svm import SVC
from sklearn.learning_curve import validation_curve

digits = load_digits()
X, y = digits.data, digits.target

param_range = np.logspace(-6, -1, 5)
train_scores, test_scores = validation_curve(
    SVC(), X, y, param_name="gamma", param_range=param_range,
    cv=10, scoring="accuracy", n_jobs=1)

plt.title("Validation Curve with SVM")
plt.xlabel("$\gamma$")
plt.ylabel("Score")
plt.ylim(0.0, 1.1)
plt.semilogx(param_range, train_scores, label="Training score")
plt.semilogx(param_range, test_scores, label="Cross-validation score")
plt.legend(loc="best")
plt.show()
