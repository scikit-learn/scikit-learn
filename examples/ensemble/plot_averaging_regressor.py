"""
===================================================
Plot individual and averaged regression predictions
===================================================

Plot individual and averaged regression predictions for Boston dataset.

First, three exemplary regressors are initialized (`GradientBoostingRegressor`,
`RandomForestRegressor`, and `LinearRegression`) and used to initialize a
`VotingRegressor`.

On plot red starred dots is averaged predictions.

"""
print(__doc__)

import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import VotingRegressor

# Loading some example data
boston = datasets.load_boston()
X = boston.data
y = boston.target

# Training classifiers
reg1 = GradientBoostingRegressor(random_state=1, n_estimators=10)
reg2 = RandomForestRegressor(random_state=1, n_estimators=10)
reg3 = LinearRegression()
ereg = VotingRegressor([('gb', reg1), ('rf', reg2), ('lr', reg3)])
reg1.fit(X, y)
reg2.fit(X, y)
reg3.fit(X, y)
ereg.fit(X, y)

xt = X[:20]
yt = y[:20]

plt.figure()
plt.plot(yt, reg1.predict(xt), 'gd', label='GradientBoostingRegressor')
plt.plot(yt, reg2.predict(xt), 'b^', label='RandomForestRegressor')
plt.plot(yt, reg3.predict(xt), 'ys', label='LinearRegression')
plt.plot(yt, ereg.predict(xt), 'r*', label='VotingRegressor')
plt.ylabel('predicted')
plt.xlabel('experimental')
plt.legend(loc="best")
plt.show()
