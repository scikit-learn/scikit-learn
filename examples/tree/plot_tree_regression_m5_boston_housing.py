"""
A copy of the Scikit Learn Gradient Boosting regression example, with the M5P in addition

http://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_regression.html
:return:
"""
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.tree import M5Prime
from sklearn.utils import shuffle
from sklearn import datasets, ensemble

import matplotlib.pyplot as plt
# plt.ion()  # comment to debug

# #############################################################################
# Load data
boston = datasets.load_boston()
X, y = shuffle(boston.data, boston.target, random_state=13)
X = X.astype(np.float32)
offset = int(X.shape[0] * 0.9)
X_train, y_train = X[:offset], y[:offset]
X_test, y_test = X[offset:], y[offset:]

# #############################################################################
# Fit regression model
params = {'n_estimators': 500, 'max_depth': 4, 'min_samples_split': 2,
          'learning_rate': 0.01, 'loss': 'ls'}
clf = ensemble.GradientBoostingRegressor(**params)
clf.fit(X_train, y_train)
y_predicted = clf.predict(X_test)
mse = mean_squared_error(y_test, y_predicted)
print("XGBoost MSE: %.4f" % mse)

# #############################################################################
# Plot predictions
plt.figure(figsize=(18, 12))
plt.subplot(2, 3, 1)
plt.title('Predictions on test set (RMSE = {:2f})'.format(np.sqrt(mse)))
plt.plot(y_test, y_predicted, '.')
plt.xlabel('true y')
plt.ylabel('predicted_y')

# #############################################################################
# Plot training deviance

# compute test set deviance
test_score = np.zeros((params['n_estimators'],), dtype=np.float64)

for i, y_pred in enumerate(clf.staged_predict(X_test)):
    test_score[i] = clf.loss_(y_test, y_pred)

plt.subplot(2, 3, 2)
plt.title('Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, clf.train_score_, 'b-',
         label='Training Set Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
         label='Test Set Deviance')
plt.legend(loc='upper right')
plt.xlabel('Boosting Iterations')
plt.ylabel('Deviance')

# #############################################################################
# Plot feature importance
feature_importance = clf.feature_importances_
# make importances relative to max importance
feature_importance = 100.0 * (feature_importance / feature_importance.max())
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
plt.subplot(2, 3, 3)
plt.barh(pos, feature_importance[sorted_idx], align='center')
plt.yticks(pos, boston.feature_names[sorted_idx])
plt.xlabel('Relative Importance')
plt.title('Variable Importance')

# ------- M5P
# #############################################################################
# Fit regression model
params = {}
clf = M5Prime(**params)
clf.fit(X_train, y_train)

# Print the tree TODO
# print(tree_to_text(clf, out_file=None))
# print(tree_to_text(clf, out_file=None, feature_names=boston.feature_names))

# Predict
y_predicted = clf.predict(X_test)
mse = mean_squared_error(y_test, y_predicted)
print("M5P MSE: %.4f" % mse)

# #############################################################################
# Plot predictions
plt.subplot(2, 3, 4)
plt.title('Predictions on test set (RMSE = {:2f})'.format(np.sqrt(mse)))
plt.plot(y_test, y_predicted, '.')
plt.xlabel('true y')
plt.ylabel('predicted_y')

# Compress the tree (features-wise)
idx = clf.compress_features()
# Print the tree TODO
# print(tree_to_text(clf, out_file=None))
# print(tree_to_text(clf, out_file=None, feature_names=boston.feature_names))

# Predict
y_predicted2 = clf.predict(X_test[:, idx])
mse2 = mean_squared_error(y_test, y_predicted2)
print("M5P2 MSE: %.4f" % mse)

# #############################################################################
# Plot predictions
plt.subplot(2, 3, 5)
plt.title('Predictions on test set (RMSE = {:2f})'.format(np.sqrt(mse2)))
plt.plot(y_test, y_predicted2, '.')
plt.xlabel('true y')
plt.ylabel('predicted_y')

# #############################################################################
# Plot feature importance
feature_importance = clf.feature_importances_
# make importances relative to max importance
feature_importance = 100.0 * (feature_importance / feature_importance.max())
sorted_idx = np.argsort(feature_importance)
pos = np.arange(sorted_idx.shape[0]) + .5
plt.subplot(2, 3, 6)
plt.barh(pos, feature_importance[sorted_idx], align='center')
# do not forget that we now work on reindexed features
plt.yticks(pos, boston.feature_names[idx][sorted_idx])
plt.xlabel('Relative Importance')
plt.title('Variable Importance')

plt.show()
