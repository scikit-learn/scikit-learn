"""
=================================================
Plot individual and voting regression predictions
=================================================

.. currentmodule:: sklearn

A voting regressor is an ensemble meta-estimator that fits base regressors each
on the whole dataset. It, then, averages the individual predictions to form a
final prediction.
We will use three different regressors to predict the data:
:class:`~ensemble.GradientBoostingRegressor`,
:class:`~ensemble.RandomForestRegressor`, and
:class:`~linear_model.LinearRegression`).
Then, using them we will make voting regressor
:class:`~ensemble.VotingRegressor`.

Finally, we will plot all of them for comparison.

We will work with the diabetes dataset which consists of the 10 features
collected from a cohort of diabetes patients. The target is the disease
progression after one year from the baseline.

"""
print(__doc__)

import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import VotingRegressor

##############################################################################
# Training classifiers
# --------------------------------
#
# First, we are going to load diabetes dataset and initiate gradient boosting
# regressor, random forest regressor and linear regression. Next, we are going
# to use each of them to build the voting regressor:

X, y = datasets.load_diabetes(return_X_y=True)

# Train classifiers
reg1 = GradientBoostingRegressor(random_state=1)
reg2 = RandomForestRegressor(random_state=1)
reg3 = LinearRegression()

reg1.fit(X, y)
reg2.fit(X, y)
reg3.fit(X, y)

ereg = VotingRegressor([('gb', reg1), ('rf', reg2), ('lr', reg3)])
ereg.fit(X, y)

##############################################################################
# Making predictions
# --------------------------------
#
# Now we will use each of the regressors to make 20 first predictions about the
# diabetes dataset.

xt = X[:20]

pred1 = reg1.predict(xt)
pred2 = reg2.predict(xt)
pred3 = reg3.predict(xt)
pred4 = ereg.predict(xt)

##############################################################################
# Plot the results
# --------------------------------
#
# Finally, we will visualize the 20 predictions. The red stars show the average
# prediction

plt.figure()
plt.plot(pred1, 'gd', label='GradientBoostingRegressor')
plt.plot(pred2, 'b^', label='RandomForestRegressor')
plt.plot(pred3, 'ys', label='LinearRegression')
plt.plot(pred4, 'r*', ms=10, label='VotingRegressor')

plt.tick_params(axis='x', which='both', bottom=False, top=False,
                labelbottom=False)
plt.ylabel('predicted')
plt.xlabel('training samples')
plt.legend(loc="best")
plt.title('Regressor predictions and their average')

plt.show()
