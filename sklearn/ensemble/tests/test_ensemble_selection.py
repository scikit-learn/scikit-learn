"""Test the ensemble_selection classifier and regressor."""

import numpy as np
from sklearn.datasets import load_iris
from sklearn.datasets import load_diabetes
from sklearn.model_selection import train_test_split


from sklearn.svm import LinearSVR
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neighbors import KNeighborsRegressor
from sklearn.neural_network import MLPClassifier
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process import GaussianProcessRegressor

from sklearn import preprocessing
from sklearn.preprocessing import scale
from sklearn.metrics import log_loss, mean_squared_error

from sklearn.ensemble import EnsembleSelection


# Regression class
diabetes = load_diabetes()
X_diabetes, y_diabetes = diabetes.data, diabetes.target

# Classification task
iris = load_iris()
X_iris, y_iris = iris.data, iris.target


def test_ensemble_selection_regressor():
    seed = 1
    X_train, X_test, y_train, y_test = train_test_split(
        scale(X_diabetes), y_diabetes, random_state=seed
    )

    # build estimators
    estimators = []
    estimators.append(("est1", GaussianProcessRegressor(random_state=seed)))
    estimators.append(("est2", RandomForestRegressor(random_state=seed)))
    estimators.append(("est3", KNeighborsRegressor()))
    estimators.append(("est4", MLPClassifier(random_state=seed)))
    estimators.append(("est5", LinearSVR(random_state=seed)))

    # fit estimators
    indiv_score = []
    for _, base_model in estimators:
        base_model.fit(X_train, y_train)
        y_pred = base_model.predict(X_test)
        loss = mean_squared_error(y_test, y_pred)
        indiv_score.append(loss)

    # fit ensemble_selector
    clf = EnsembleSelection(
        estimators=estimators,
        score_metric=mean_squared_error,
        score_direction="min",
        is_base_estimator_proba=False,
        min_estimators=2,
        max_estimators=10,
        with_replacement=True,
        verbose=True,
    )
    clf.fit(X_test, y_test)
    y_pred = clf.predict(X_test)

    loss = mean_squared_error(y_test, y_pred)
    assert loss <= np.min(indiv_score)


def test_ensemble_selection_classifier():
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    seed = 3

    X_train, X_test, y_train, y_test = train_test_split(
        scale(X_iris), y_iris, stratify=y_iris, test_size=0.5, random_state=seed
    )

    # convert label id to one-hot vector
    y_test_one_hot = (
        preprocessing.OneHotEncoder().fit_transform(y_test.reshape((-1, 1))).toarray()
    )

    # instantiate estimators
    estimators = []
    estimators.append(("est1", GaussianProcessClassifier(random_state=seed)))
    estimators.append(("est2", SVC(random_state=seed, probability=True)))
    estimators.append(("est3", RandomForestClassifier(random_state=seed)))
    estimators.append(("est4", KNeighborsClassifier()))
    estimators.append(("est5", MLPClassifier(random_state=seed)))

    # fit estimators
    for name, base_model in estimators:
        base_model.fit(X_train, y_train)

    clf = EnsembleSelection(
        estimators=estimators,
        score_metric=log_loss,
        score_direction="min",
        is_base_estimator_proba=True,
        min_estimators=2,
        max_estimators=10,
        with_replacement=True,
        verbose=True,
    )
    clf.fit(X_test, y_test_one_hot)
    y_pred = clf.predict_proba(X_test)

    loss = log_loss(y_test, y_pred)
    assert loss < 0.2  # not so bad


def test_ensemble_selection_classifier_unamed_estimators():
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    seed = 3

    X_train, X_test, y_train, y_test = train_test_split(
        scale(X_iris), y_iris, stratify=y_iris, test_size=0.5, random_state=seed
    )

    # convert label id to one-hot vector
    y_test_one_hot = (
        preprocessing.OneHotEncoder().fit_transform(y_test.reshape((-1, 1))).toarray()
    )

    # instantiate estimators
    estimators = []
    estimators.append((GaussianProcessClassifier(random_state=seed)))
    estimators.append((SVC(random_state=seed, probability=True)))
    estimators.append((RandomForestClassifier(random_state=seed)))
    estimators.append((KNeighborsClassifier()))
    estimators.append((MLPClassifier(random_state=seed)))

    # fit estimators
    for base_model in estimators:
        base_model.fit(X_train, y_train)

    clf = EnsembleSelection(
        estimators=estimators,
        score_metric=log_loss,
        score_direction="min",
        is_base_estimator_proba=True,
        min_estimators=2,
        max_estimators=10,
        with_replacement=True,
        verbose=True,
    )
    clf.fit(X_test, y_test_one_hot)
    y_pred = clf.predict_proba(X_test)

    loss = log_loss(y_test, y_pred)
    assert loss < 0.2  # not so bad


def test_ensemble_selection_classifier_sparse_labels_given():
    # prescale the data to avoid convergence warning without using a pipeline
    # for later assert
    seed = 3

    X_train, X_test, y_train, y_test = train_test_split(
        scale(X_iris), y_iris, stratify=y_iris, test_size=0.5, random_state=seed
    )

    # instantiate estimators
    estimators = []
    estimators.append(("A", GaussianProcessClassifier(random_state=seed)))
    estimators.append(("B", SVC(random_state=seed, probability=True)))
    estimators.append(("C", RandomForestClassifier(random_state=seed)))
    estimators.append(("D", KNeighborsClassifier()))
    estimators.append(("E", MLPClassifier(random_state=seed)))

    # fit estimators
    for _, base_model in estimators:
        base_model.fit(X_train, y_train)

    clf = EnsembleSelection(
        estimators=estimators,
        score_metric=log_loss,
        score_direction="min",
        is_base_estimator_proba=True,
        min_estimators=2,
        max_estimators=10,
        with_replacement=True,
        verbose=True,
    )
    clf.fit(X_test, y_test)
    y_pred = clf.predict_proba(X_test)

    loss = log_loss(y_test, y_pred)
    assert loss < 0.2  # not so bad
