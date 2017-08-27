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

MSE for lasso: 22.91279
MSE for ridge: 21.69905
MSE for svr: 32.25357
MSE for Stacked Regressors: 18.93744
"""

# utils
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

# setup plot
SUBPLOT_X = 1
SUBPLOT_Y = 4
SUBPLOT_OFFSET = (SUBPLOT_X * 10 + SUBPLOT_Y) * 10
EXAMPLE_INDEX = 1

plt.figure(EXAMPLE_INDEX, figsize=(9, 3))

# prepare data
X, y = load_boston(return_X_y=True)
X_train, X_test, y_train, y_test = train_test_split(X, y,
                                                    random_state=RANDOM_SEED)

# Base regressors
lasso = LassoCV(random_state=RANDOM_SEED)
ridge = RidgeCV()
svr = SVR(C=1e2, gamma=1e-3)

base_regressors = [("lasso", lasso),
                   ("ridge", ridge),
                   ("svr", svr)]


def evaluate_and_log_model(name, model, plot_idx):
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    score = mean_squared_error(y_test, y_pred)
    print("MSE for %s: %.5f" % (name, score))

    plt.subplot(SUBPLOT_OFFSET + plot_idx)
    plt.title(name)
    plt.scatter(y_pred, y_test, )
    plt.xlabel('Predicted')
    plt.ylabel('Measured')
    plt.xticks(())
    plt.yticks(())


for i, (name, regressor) in enumerate(base_regressors):
    evaluate_and_log_model(name, regressor, 1 + i)


# Stacked ensemble: we use the base regressors as features for a new linear
# regressor.
layer0 = make_stack_layer(base_regressors)
final_regressor = Pipeline([('layer0', layer0),
                            ('layer1', LinearRegression())])

evaluate_and_log_model('Stacked Regressors', final_regressor, 4)
plt.tight_layout()
plt.show()
