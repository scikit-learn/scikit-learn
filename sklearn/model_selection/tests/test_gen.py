import pytest
import numpy as np
from sklearn.svm import SVC

from sklearn import datasets
from sklearn.model_selection import cross_validate, GeneticSearchCV


def test_gen():
    iris = datasets.load_iris()
    model = SVC()
    E = np.random.uniform(0, 0.1, size=(len(iris.data), 20))  # add some noise
    X = np.hstack((iris.data, E))
    y = iris.target
    params = dict(kernel=['linear', 'poly', 'rbf', 'sigmoid'], C=dict(min_=0.0, max_=0.1),
                  gamma=dict(min_=0.0, max_=0.1))
    cv = GeneticSearchCV(model, params, return_train_score=True, cv=3)
    result = cv.fit(X, y)
    # print(result)
    print(cv.best_score_, cv.best_params_)
