"""
=================================
Combine predictors using stacking
=================================

Stacking refers to a method to blend estimators. In this strategy, some
estimators are individually fitted on some training data while a final
estimator is trained using the stacked predictions of these base estimators.

In this example, we illustrate the use case in which different regressors
(linear and non-linear) are stacked together and a final non-linear regressor
is used to output the prediction. We compare the performance of each individual
regressor with the stacking strategy. Stacking improves the overall
performance.

"""

# Authors: Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

import numpy as np

###############################################################################
# Stack of predictors on a single data set
###############################################################################
# Linear (ridge and lasso) and non-linear (RBF-SVR) are stacked together using
# a gradient-boosting regressor.

from sklearn.linear_model import LassoCV, RidgeCV
from sklearn.svm import SVR
from sklearn.ensemble import StackingRegressor, GradientBoostingRegressor

estimators = [('Lasso', LassoCV(random_state=42)),
              ('Ridge', RidgeCV()),
              ('SVR', SVR(C=1, gamma=1e-6, kernel='rbf'))]
stacking_regressor = StackingRegressor(
    estimators=estimators,
    final_estimator=GradientBoostingRegressor(random_state=0))

###############################################################################
# The function ``plot_regression_results`` is used to plot the predicted and
# true targets.

import matplotlib.pyplot as plt


def plot_regression_results(ax, y_true, y_pred, title, scores):
    """Scatter plot of the predicted vs true targets."""
    ax.plot([y_true.min(), y_true.max()],
            [y_true.min(), y_true.max()],
            '--r', linewidth=2)
    ax.scatter(y_true, y_pred, alpha=0.2)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['bottom'].set_position(('outward', 10))
    ax.set_xlim([y_true.min(), y_true.max()])
    ax.set_ylim([y_true.min(), y_true.max()])
    ax.set_xlabel('Measured')
    ax.set_ylabel('Predicted')
    extra = plt.Rectangle((0, 0), 0, 0, fc="w", fill=False,
                          edgecolor='none', linewidth=0)
    ax.legend([extra], [scores], loc='upper left')
    ax.set_title(title)


###############################################################################
# We used the Boston data set (prediction of house prices). We check the
# performance of each individual predictor as well as the stack of the
# regressors.

from sklearn.datasets import load_boston
from sklearn.model_selection import cross_validate, cross_val_predict

X, y = load_boston(return_X_y=True)

fig, axs = plt.subplots(2, 2, figsize=(9, 7))
axs = np.ravel(axs)

for ax, (name, est) in zip(axs, estimators + [('Stacking Regressor',
                                               stacking_regressor)]):
    score = cross_validate(est, X, y,
                           scoring=['r2', 'neg_mean_absolute_error'],
                           n_jobs=-1, verbose=0)

    y_pred = cross_val_predict(est, X, y, n_jobs=-1, verbose=0)
    plot_regression_results(
        ax, y, y_pred,
        name,
        r'$R^2={:.2f} \pm {:.2f}$\n$MAE={:.2f} \pm {:.2f}$'
        .format(np.mean(score['test_r2']),
                np.std(score['test_r2']),
                np.mean(score['test_neg_mean_absolute_error']),
                np.std(score['test_neg_mean_absolute_error'])))

plt.suptitle('Single predictors versus stacked predictors')
plt.tight_layout()
plt.subplots_adjust(top=0.9)
plt.show()
