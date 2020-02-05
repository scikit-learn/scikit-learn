"""
=================================================
Plot individual and voting regression predictions
=================================================

.. currentmodule:: sklearn

Regression is a method to predict the average value of the target y from the
given data X. It is possible to use different type of regressors to predict the
data. In this example we are going to use three of them:
:class:`~ensemble.GradientBoostingRegressor`,
:class:`~ensemble.RandomForestRegressor`, and
:class:`~linear_model.LinearRegression`).
Then, we are going to use their predictions to make another one: voting
regression predictions using
:class:`~ensemble.VotingRegressor`.

Finally, we will plot all of them for comparison.

We will work with the diabetes dataset which consists of the 10 features
collected from the diabetes patients. The target is the disease progression
after one year from the baseline.

"""
print(__doc__)

import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import VotingRegressor

"""
=================================
Training classifiers
=================================
First, we are going to load diabetes dataset and initiate gradient boosting
regressor, random forest regressor and linear regression. Next, we are going to
use each of the initialized regressors to build the voting regressor.

"""

# Loading some example data
X, y = datasets.load_diabetes(return_X_y=True)

# Training classifiers
reg1 = GradientBoostingRegressor(random_state=1)
reg2 = RandomForestRegressor(random_state=1)
reg3 = LinearRegression()

reg1.fit(X, y)
reg2.fit(X, y)
reg3.fit(X, y)

ereg = VotingRegressor([('gb', reg1), ('rf', reg2), ('lr', reg3)])
ereg.fit(X, y)

"""
=================================
Making predictions
=================================
Now we will use each of the regressors to make predictions about the diabetes
dataset.

"""

pred1 = reg1.predict(X)
pred2 = reg2.predict(X)
pred3 = reg3.predict(X)
pred4 = ereg.predict(X)

"""
=================================
Plot figure
=================================
Finally, we will calculate each prediction, plot them and visualize first 20
predictions.

"""

plt.figure()
plt.plot(pred1, 'gd', label='GradientBoostingRegressor')
plt.plot(pred2, 'b^', label='RandomForestRegressor')
plt.plot(pred3, 'ys', label='LinearRegression')
plt.plot(pred4, 'r*', ms=10, label='VotingRegressor')

plt.tick_params(axis='x', which='both', bottom=False, top=False,
                labelbottom=False)
plt.ylabel('predicted')
plt.xlabel('training samples')
plt.xlim([-0.5, 20.5])
plt.legend(loc="best")
plt.title('Regressor predictions and their average')

plt.show()
