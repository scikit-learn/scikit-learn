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
:class:`~linear_model.LinearRegression`). Next, we are going to use their
predictions to make andother one: voting regression predictions using
:class:`~ensemble.VotingRegressor`.

Finally, we will plot all of them for comparison:

We will work with the diabetes dataset which consists of the 10 featuers
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
regressor, random forest regressor and linear regression. We set random_state
to 1 to have the same results after each run. Next, we are going to use each of
the initialized regressors to build the voting regressor.

"""

# Loading some example data
X, y = datasets.load_diabetes(return_X_y=True)

# Training classifiers
reg1 = GradientBoostingRegressor(random_state=1)
reg2 = RandomForestRegressor(random_state=1)
reg3 = LinearRegression()
ereg = VotingRegressor([('gb', reg1), ('rf', reg2), ('lr', reg3)])
reg1.fit(X, y)
reg2.fit(X, y)
reg3.fit(X, y)
ereg.fit(X, y)

xt = X[:20]

"""
=================================
Plot figure
=================================



"""

plt.figure()
plt.plot(reg1.predict(xt), 'gd', label='GradientBoostingRegressor')
plt.plot(reg2.predict(xt), 'b^', label='RandomForestRegressor')
plt.plot(reg3.predict(xt), 'ys', label='LinearRegression')
plt.plot(ereg.predict(xt), 'r*', label='VotingRegressor')
plt.tick_params(axis='x', which='both', bottom=False, top=False,
                labelbottom=False)
plt.ylabel('predicted')
plt.xlabel('training samples')
plt.legend(loc="best")
plt.title('Comparison of individual predictions with averaged')
plt.show()
