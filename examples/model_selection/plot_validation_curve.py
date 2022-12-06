"""
==========================
Plotting Validation Curves
==========================

In this plot you can see the training scores and validation scores of an SVM
for different values of the kernel parameter gamma. For very low values of
gamma, you can see that both the training score and the validation score are
low. This is called underfitting. Medium values of gamma will result in high
values for both scores, i.e. the classifier is performing fairly well. If gamma
is too high, the classifier will overfit, which means that the training score
is good but the validation score is poor.

"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import load_digits
from sklearn.svm import SVC
from sklearn.model_selection import ValidationCurveDisplay

X, y = load_digits(return_X_y=True)
subset_mask = np.isin(y, [1, 2])  # binary classification: 1 vs 2
X, y = X[subset_mask], y[subset_mask]

param_name, param_range = "gamma", np.logspace(-6, -1, 5)
disp = ValidationCurveDisplay.from_estimator(
    SVC(),
    X,
    y,
    param_name=param_name,
    param_range=param_range,
    score_type="both",
    log_scale=True,
    n_jobs=2,
    score_name="Accuracy",
)
disp.ax_.set_title("Validation Curve with SVM")
disp.ax_.set_xlabel(r"$\gamma$")
disp.ax_.set_ylim(0.0, 1.1)
plt.show()
