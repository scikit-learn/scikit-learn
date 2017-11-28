"""
======================================================
Stacked generalization applied to a regression problem
======================================================

This example shows how stacked generalization can be used to combine several
regressors into a single stacked one that performs better than the best
regressor.

We use Boston's house pricing dataset to compare the mean squared error between
three regressors (SVM, Lasso and Ridge regressions) and the combination of
their outputs with a single linear regression. The following result is
achieved.
"""

# utils
from time import time
import matplotlib.pyplot as plt

# import base regressors
from sklearn.linear_model import LassoCV, RidgeCV, LinearRegression
from sklearn.svm import SVR

# stacking api
from sklearn.ensemble import make_stack_layer
from sklearn.pipeline import Pipeline

# dataset
from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error

RANDOM_SEED = 89555

fig, axarr = plt.subplots(1, 4, figsize=(9, 3))

# prepare data
X, y = load_boston(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state=RANDOM_SEED)

# Base regressors
lasso = LassoCV(random_state=RANDOM_SEED)
ridge = RidgeCV()
svr = SVR(C=1e3, gamma=1e-4, kernel='rbf')

base_regressors = [("Lasso Regressor", lasso),
                   ("Ridge Regressor", ridge),
                   ("SVR", svr)]


def evaluate_and_log_model(name, model, ax):
    t0_train = time()
    model.fit(X_train, y_train)
    train_time = time() - t0_train
    y_pred = model.predict(X_test)
    score = mean_squared_error(y_test, y_pred)
    print("MSE for %s: %.3f (train time: %.3f seconds)"
          % (name, score, train_time))

    ax.scatter(y_test, y_pred, edgecolors=(0, 0, 0))
    ax.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()],
            'k--', lw=4)
    ax.set_title(name)
    ax.set_xlabel('Measured')
    ax.set_ylabel('Predicted')
    ax.set_xticks(())
    ax.set_yticks(())


for i, (name, regressor) in enumerate(base_regressors):
    evaluate_and_log_model(name, regressor, axarr[i])

# Stacked ensemble: we use the base regressors as features for a new linear
# regressor.
layer0 = make_stack_layer(base_regressors)
final_regressor = Pipeline([('layer0', layer0),
                            ('layer1', LinearRegression())])

evaluate_and_log_model('Stacked Regressors', final_regressor, axarr[-1])
fig.tight_layout()
fig.show()
