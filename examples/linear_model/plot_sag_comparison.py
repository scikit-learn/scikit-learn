"""
=========================================
Comparing SAG to plain SGD implementation
=========================================

An example comparing SAG with SGD

"""
# Author: Danny Sullivan <dsullivan7 at hotmail dot com>

import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets

from sklearn.cross_validation import train_test_split
from sklearn.linear_model import (SGDClassifier, SAGClassifier,
                                  SGDRegressor, SAGRegressor)
from sklearn.metrics import mean_squared_error

heldout = [0.95, 0.90, 0.75, 0.50, 0.25, 0.10, 0.01]
rounds = 5
digits = datasets.load_digits()
digits.data = digits.data[digits.target < 4]
digits.target = digits.target[digits.target < 4]
boston = datasets.load_boston()

alpha = .0001
classifiers = [
    ("SGD Classifier", SGDClassifier(loss="log", alpha=alpha)),
    ("SAG Classifier", SAGClassifier(eta0='auto', alpha=alpha,
                                     random_state=77)),
]

alpha = .0001
regressors = [
    ("SGD Regressor", SGDRegressor(eta0=.00001, loss="squared_loss",
                                   alpha=alpha)),
    ("SAG Regressor", SAGRegressor(eta0='auto', alpha=alpha, random_state=77)),
]

xx = 1. - np.array(heldout)


def c_score(clf, X, y):
    return 1 - clf.score(X, y)


def r_score(clf, X, y):
    return mean_squared_error(y, clf.predict(X))

all_tests = (
    (boston, regressors, r_score, "Mean Squared Error"),
    (digits, classifiers, c_score, "Accuracy Error")
)

plt.figure()
for data_set, clfs, scoring_func, plot_label in all_tests:
    X = data_set.data
    y = data_set.target
    for name, clf in clfs:
        print("training", name, "...")
        rng = np.random.RandomState(42)
        yy = []
        for i in heldout:
            yy_ = []
            for r in range(rounds):
                X_train, X_test, y_train, y_test = \
                    train_test_split(X, y, test_size=i, random_state=rng)
                clf.fit(X_train, y_train)
                yy_.append(scoring_func(clf, X_test, y_test))
            yy.append(np.mean(yy_))
        plt.plot(xx, yy, label=name)

    plt.legend(loc="upper right")
    plt.xlabel("Proportion train")
    plt.ylabel(plot_label)
    plt.show()
