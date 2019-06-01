"""
==============================
Dummy baselines usage examples
==============================

Example of DummyClassifier and DummyRegression.

A model or pipeline's result will typically be compared to a baseline. Simple
baselines are provided for both classification and regression problems through
the DummyClassifier and DummyRegression estimators.
Therefore it is expected that a good model or pipeline should beat the dummy
heuristics implemented here. The only way to find it out is to run those
dummy estimators prior to modelling your pipeline.

"""
print(__doc__)

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.dummy import DummyRegressor, DummyClassifier
from sklearn.linear_model import LinearRegression, LogisticRegression
from sklearn.datasets import make_regression

print("*" * 80)
print("Regression example")
print("*" * 80)
print("Since all features are informative the LinearRegression model should")
print("capture perfect score while the DummyRegressor should give a near")
print("non-X dependent result (a poor 0 score)")


X, y = make_regression(random_state=0, n_samples=1000, n_features=100,
                       n_informative=100)

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

clf = DummyRegressor(strategy="mean")
clf.fit(X_train, y_train)
print("Dummy mean score {}".format(clf.score(X_test, y_test)))

clf = LinearRegression()
clf.fit(X_train, y_train)
print("Linear regression score {}".format(clf.score(X_test, y_test)))


print("*" * 80)
print("Time series regression example")
print("*" * 80)
print("Since this is a fibonnaci alike series, the last element is a good")
print("predictor of the next one, giving the dummy regressor an informative")
print("score (in this case a score between 0 and 1).")
print("But the linear regression should capture the full relationship")
print("between N and its two previous elements (score near 1).")

# sample time series with a linear relationship N = (N-1) + (N-2)
full_time_series = [1, 2, 3, 5, 8, 13, 21, 34, 55, 87, 142]

# splits the time series in a rolling window of size 2
steps_for_prediction = 2
count = len(full_time_series) - steps_for_prediction
X = [full_time_series[i:i+steps_for_prediction] for i in range(count)]

# y is the next step
y = full_time_series[steps_for_prediction:]

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)

clf = DummyRegressor(strategy="series_last")
clf.fit(X_train, y_train)
print("Dummy score {}".format(clf.score(X_test, y_test)))

clf = LinearRegression()
clf.fit(X_train, y_train)
print("Linear regression score {}".format(clf.score(X_test, y_test)))


print("*" * 80)
print("Classification example")
print("*" * 80)
print("Since the relationship between y and X is linear ")
print("(X > 500 => 1), logistic regression should give a")
print("near perfect score (1) while the default")
print("DummyClassifier will give a near 0.5 result.")

X = [[x] for x in range(1000)]
y = np.array([x > 500 for x in range(1000)]) * 1

X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0,
                                                    stratify=y)

clf = DummyClassifier(random_state=0)
clf.fit(X_train, y_train)
print("Dummy score {}".format(clf.score(X_test, y_test)))

clf = LogisticRegression(random_state=0)
clf.fit(X_train, y_train)
print("Logistic regression score {}".format(clf.score(X_test, y_test)))
