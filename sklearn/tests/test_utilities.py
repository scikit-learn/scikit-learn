# -*- coding: utf8 -*-
# Author: Alain Pena <alain.pena.r@gmail.com>
#
# License: BSD 3 clause

import numpy as np
from itertools import izip_longest

from sklearn.datasets import make_multilabel_classification
from sklearn.utils.validation import check_random_state
from sklearn.utils import testing as testing
from sklearn.multiclass import OneVsRestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn import utilities as util


def test_convert_probabilities():
    n_features = 10
    n_labels = 2
    n_samples = 10
    n_classes = 5
    length = 1
    random_state = 0
    X, y = make_multilabel_classification(
        n_samples=n_samples, n_features=n_features, length=length,
        n_classes=n_classes, n_labels=n_labels,
        allow_unlabeled=False,
        random_state=random_state, return_indicator=True)
    Xt, yt = make_multilabel_classification(
        n_samples=n_samples, n_features=n_features, length=length,
        n_classes=n_classes, n_labels=n_labels,
        allow_unlabeled=False,
        random_state=random_state + n_samples, return_indicator=True)
    tree = DecisionTreeClassifier(random_state=random_state)
    tree.fit(X, y)
    base_mm = tree.predict_proba(Xt)
    base_1m = tree.predict_proba(Xt[0, :].reshape(1, -1))
    tree.fit(X, y[:, 0].ravel())
    base_m1 = tree.predict_proba(Xt)
    base_11 = tree.predict_proba(Xt[0, :].reshape(1, -1))

    trees = ExtraTreesClassifier(random_state=random_state)
    onevs = OneVsRestClassifier(tree)
    onevss = OneVsRestClassifier(trees)

    def do_one(estimator):
        print(estimator)
        estimator.fit(X, y)
        x = util.convert_probabilities(estimator.predict_proba(Xt), y.shape)
        testing.assert_equal(type(base_mm), type(x))
        testing.assert_equal(len(base_mm), len(x))
        for a, b in zip(base_mm, x):
            testing.assert_array_equal(a.shape, b.shape)
        c = util.convert_probabilities(x, y.shape)
        testing.assert_equal(type(c), type(x))
        testing.assert_equal(len(c), len(x))
        for a, b in zip(c, x):
            testing.assert_array_equal(a.shape, b.shape)

        x = util.convert_probabilities(
            estimator.predict_proba(Xt[0, :].reshape(1, -1)), y.shape)
        testing.assert_equal(type(base_1m), type(x))
        testing.assert_equal(len(base_1m), len(x))
        for a, b in zip(base_1m, x):
            testing.assert_array_equal(a.shape, b.shape)
        c = util.convert_probabilities(x, y.shape)
        testing.assert_equal(type(c), type(x))
        testing.assert_equal(len(c), len(x))
        for a, b in zip(c, x):
            testing.assert_array_equal(a.shape, b.shape)

        estimator.fit(X, y[:, 0].ravel())
        x = util.convert_probabilities(
            estimator.predict_proba(Xt), y[:, 0].ravel().shape)
        testing.assert_equal(type(base_m1), type(x))
        testing.assert_equal(base_m1.shape, x.shape)
        c = util.convert_probabilities(x, y[:, 0].ravel().shape)
        testing.assert_equal(type(c), type(x))
        testing.assert_equal(len(c), len(x))
        for a, b in zip(c, x):
            testing.assert_array_equal(a.shape, b.shape)

        x = util.convert_probabilities(
            estimator.predict_proba(Xt[0, :].reshape(1, -1)),
            y[:, 0].ravel().shape)
        testing.assert_equal(type(base_11), type(x))
        testing.assert_equal(base_11.shape, x.shape)
        c = util.convert_probabilities(x, y[:, 0].ravel().shape)
        testing.assert_equal(type(c), type(x))
        testing.assert_equal(len(c), len(x))
        for a, b in zip(c, x):
            testing.assert_array_equal(a.shape, b.shape)
    do_one(trees)
    do_one(onevs)
    do_one(onevss)


def test_cycle_permutations():
    l = []
    res = [['2', '1', '0'],
           ['0', '1', '2'],
           ['0', '2', '1'],
           ['1', '0', '2'],
           ['1', '2', '0'],
           ['2', '0', '1'],
           ['2', '1', '0'],
           ['0', '1', '2'],
           ['0', '2', '1'],
           ['1', '0', '2']]
    for p in util.cycle_permutations("201", n_iter=10):
        l.append(p)
    testing.assert_array_equal(l, res)
    l = []
    for p in util.cycle_permutations("210", n_iter=9):
        l.append(p)
    testing.assert_array_equal(l, res[1:])
    l = []
    for p in util.cycle_permutations("210", n_iter=0):
        l.append(p)
    testing.assert_array_equal(l, [])
    for p in util.cycle_permutations("210", n_iter=-5):
        l.append(p)
    testing.assert_array_equal(l, [])
    l = []
    for p in util.cycle_permutations("201", n_iter=1):
        l.append(p)
    testing.assert_array_equal(l, [['2', '1', '0']])
    l = 0
    for p in util.cycle_permutations("210", n_iter=None):
        l = l + 1
        if (l > 100):
            break
    l = []
    for p in util.cycle_permutations(iter(res[0]), n_iter=9):
        l.append(p)
    testing.assert_array_equal(l, res[1:])
    for p in util.cycle_permutations(32541, n_iter=10):
        raise AssertionError("Error: cycle_permutations "
                             "failed for int(32541).")
    for idx, p in enumerate(util.cycle_permutations("0", n_iter=10)):
        testing.assert_array_equal(p, ['0'])
        if (idx > 8):
            break
    else:
        raise AssertionError("Error: cycle_permutations failed for '0'.")
    for idx, p in enumerate(util.cycle_permutations("000", n_iter=10)):
        testing.assert_array_equal(p, ['0', '0', '0'])
        if (idx > 8):
            break
    else:
        raise AssertionError("Error: cycle_permutations failed for '000'.")


def test_get_most_probable_class():
    test = [[0.3, 0.07, 0.63],
            [1.0, 0.00, 0.00],
            [0.4, 0.60, 0.00],
            [0.5, 0.50, 0.00]]
    testing.assert_array_equal(util.get_most_probable_class(test),
                               [2, 0, 1, 0])
    test = [[0.1, 0.2, 0.5, 0.2],
            [0.4, 0.1, 0.1, 0.4],
            [0.1, 0.4, 0.2, 0.3]]
    testing.assert_array_equal(util.get_most_probable_class(test),
                               [2, 0, 1])


def test_get_random_numbers():
    shape = None
    random_state = 0
    testing.assert_equal(util.get_random_numbers(
        shape=shape, random_state=random_state), check_random_state(
        random_state).randint(util.MAX_SEED_INT, size=shape))
    testing.assert_equal(util.get_random_numbers(
        shape=shape, random_state=random_state), check_random_state(
        random_state).randint(util.MAX_SEED_INT, size=shape))
    testing.assert_not_equal(util.get_random_numbers(
        shape=shape, random_state=random_state), check_random_state(
        random_state + 1).randint(util.MAX_SEED_INT, size=shape))
    shape = (2, 2)
    testing.assert_array_equal(util.get_random_numbers(
        shape=shape, random_state=random_state), check_random_state(
        random_state).randint(util.MAX_SEED_INT, size=shape))
    testing.assert_array_equal(util.get_random_numbers(
        shape=shape, random_state=random_state), check_random_state(
        random_state).randint(util.MAX_SEED_INT, size=shape))
    try:
        testing.assert_array_equal(util.get_random_numbers(
            shape=shape, random_state=random_state), check_random_state(
            random_state + 1).randint(util.MAX_SEED_INT, size=shape))
    except AssertionError:
        pass
    else:
        testing.assert_true(
            False, msg="Error: Random generator not working for shapes.")


def test_iterize():
    testing.assert_array_equal(
        list(xrange(10)), util.iterize(list(xrange(10))))
    elem = "should not be equal"
    testing.assert_not_equal(elem, util.iterize(elem))
    for i in util.iterize(10):
        testing.assert_equal(i, 10)
    for i in util.iterize(None):
        testing.assert_equal(i, None)


def test_next_permutation():
    testing.assert_array_equal(util.next_permutation("32541"),
                               list("34125"))
    testing.assert_array_equal(util.next_permutation(("3", "2", "5",
                                                      "4", "1")),
                               list("34125"))
    testing.assert_array_equal(util.next_permutation(iter(("3", "2", "5",
                                                           "4", "1"))),
                               list("34125"))
    if not (None is util.next_permutation(32541)):
        raise AssertionError("Error: next_permutation failed for int(32541).")
    if not (None is util.next_permutation(['4', '3', '2', '2', '1', '0'])):
        raise AssertionError("Error: next_permutation failed for int(432210).")


def test_next_permutations():
    valid = [['0', '1', '2', '2', '3', '4'],
             ['0', '1', '2', '2', '4', '3'],
             ['0', '1', '2', '3', '2', '4'],
             ['0', '1', '2', '3', '4', '2'],
             ['0', '1', '2', '4', '2', '3'],
             ['0', '1', '2', '4', '3', '2'],
             ['0', '1', '3', '2', '2', '4'],
             ['0', '1', '3', '2', '4', '2'],
             ['0', '1', '3', '4', '2', '2'],
             ['0', '1', '4', '2', '2', '3'],
             ['0', '1', '4', '2', '3', '2'],
             ['0', '1', '4', '3', '2', '2'],
             ['0', '2', '1', '2', '3', '4'],
             ['0', '2', '1', '2', '4', '3'],
             ['0', '2', '1', '3', '2', '4'],
             ['0', '2', '1', '3', '4', '2'],
             ['0', '2', '1', '4', '2', '3'],
             ['0', '2', '1', '4', '3', '2'],
             ['0', '2', '2', '1', '3', '4'],
             ['0', '2', '2', '1', '4', '3'],
             ['0', '2', '2', '3', '1', '4'],
             ['0', '2', '2', '3', '4', '1'],
             ['0', '2', '2', '4', '1', '3'],
             ['0', '2', '2', '4', '3', '1'],
             ['0', '2', '3', '1', '2', '4'],
             ['0', '2', '3', '1', '4', '2'],
             ['0', '2', '3', '2', '1', '4'],
             ['0', '2', '3', '2', '4', '1'],
             ['0', '2', '3', '4', '1', '2'],
             ['0', '2', '3', '4', '2', '1'],
             ['0', '2', '4', '1', '2', '3'],
             ['0', '2', '4', '1', '3', '2'],
             ['0', '2', '4', '2', '1', '3'],
             ['0', '2', '4', '2', '3', '1'],
             ['0', '2', '4', '3', '1', '2'],
             ['0', '2', '4', '3', '2', '1'],
             ['0', '3', '1', '2', '2', '4'],
             ['0', '3', '1', '2', '4', '2'],
             ['0', '3', '1', '4', '2', '2'],
             ['0', '3', '2', '1', '2', '4'],
             ['0', '3', '2', '1', '4', '2'],
             ['0', '3', '2', '2', '1', '4'],
             ['0', '3', '2', '2', '4', '1'],
             ['0', '3', '2', '4', '1', '2'],
             ['0', '3', '2', '4', '2', '1'],
             ['0', '3', '4', '1', '2', '2'],
             ['0', '3', '4', '2', '1', '2'],
             ['0', '3', '4', '2', '2', '1'],
             ['0', '4', '1', '2', '2', '3'],
             ['0', '4', '1', '2', '3', '2'],
             ['0', '4', '1', '3', '2', '2'],
             ['0', '4', '2', '1', '2', '3'],
             ['0', '4', '2', '1', '3', '2'],
             ['0', '4', '2', '2', '1', '3'],
             ['0', '4', '2', '2', '3', '1'],
             ['0', '4', '2', '3', '1', '2'],
             ['0', '4', '2', '3', '2', '1'],
             ['0', '4', '3', '1', '2', '2'],
             ['0', '4', '3', '2', '1', '2'],
             ['0', '4', '3', '2', '2', '1'],
             ['1', '0', '2', '2', '3', '4'],
             ['1', '0', '2', '2', '4', '3'],
             ['1', '0', '2', '3', '2', '4'],
             ['1', '0', '2', '3', '4', '2'],
             ['1', '0', '2', '4', '2', '3'],
             ['1', '0', '2', '4', '3', '2'],
             ['1', '0', '3', '2', '2', '4'],
             ['1', '0', '3', '2', '4', '2'],
             ['1', '0', '3', '4', '2', '2'],
             ['1', '0', '4', '2', '2', '3'],
             ['1', '0', '4', '2', '3', '2'],
             ['1', '0', '4', '3', '2', '2'],
             ['1', '2', '0', '2', '3', '4'],
             ['1', '2', '0', '2', '4', '3'],
             ['1', '2', '0', '3', '2', '4'],
             ['1', '2', '0', '3', '4', '2'],
             ['1', '2', '0', '4', '2', '3'],
             ['1', '2', '0', '4', '3', '2'],
             ['1', '2', '2', '0', '3', '4'],
             ['1', '2', '2', '0', '4', '3'],
             ['1', '2', '2', '3', '0', '4'],
             ['1', '2', '2', '3', '4', '0'],
             ['1', '2', '2', '4', '0', '3'],
             ['1', '2', '2', '4', '3', '0'],
             ['1', '2', '3', '0', '2', '4'],
             ['1', '2', '3', '0', '4', '2'],
             ['1', '2', '3', '2', '0', '4'],
             ['1', '2', '3', '2', '4', '0'],
             ['1', '2', '3', '4', '0', '2'],
             ['1', '2', '3', '4', '2', '0'],
             ['1', '2', '4', '0', '2', '3'],
             ['1', '2', '4', '0', '3', '2'],
             ['1', '2', '4', '2', '0', '3'],
             ['1', '2', '4', '2', '3', '0'],
             ['1', '2', '4', '3', '0', '2'],
             ['1', '2', '4', '3', '2', '0'],
             ['1', '3', '0', '2', '2', '4'],
             ['1', '3', '0', '2', '4', '2'],
             ['1', '3', '0', '4', '2', '2'],
             ['1', '3', '2', '0', '2', '4'],
             ['1', '3', '2', '0', '4', '2'],
             ['1', '3', '2', '2', '0', '4'],
             ['1', '3', '2', '2', '4', '0'],
             ['1', '3', '2', '4', '0', '2'],
             ['1', '3', '2', '4', '2', '0'],
             ['1', '3', '4', '0', '2', '2'],
             ['1', '3', '4', '2', '0', '2'],
             ['1', '3', '4', '2', '2', '0'],
             ['1', '4', '0', '2', '2', '3'],
             ['1', '4', '0', '2', '3', '2'],
             ['1', '4', '0', '3', '2', '2'],
             ['1', '4', '2', '0', '2', '3'],
             ['1', '4', '2', '0', '3', '2'],
             ['1', '4', '2', '2', '0', '3'],
             ['1', '4', '2', '2', '3', '0'],
             ['1', '4', '2', '3', '0', '2'],
             ['1', '4', '2', '3', '2', '0'],
             ['1', '4', '3', '0', '2', '2'],
             ['1', '4', '3', '2', '0', '2'],
             ['1', '4', '3', '2', '2', '0'],
             ['2', '0', '1', '2', '3', '4'],
             ['2', '0', '1', '2', '4', '3'],
             ['2', '0', '1', '3', '2', '4'],
             ['2', '0', '1', '3', '4', '2'],
             ['2', '0', '1', '4', '2', '3'],
             ['2', '0', '1', '4', '3', '2'],
             ['2', '0', '2', '1', '3', '4'],
             ['2', '0', '2', '1', '4', '3'],
             ['2', '0', '2', '3', '1', '4'],
             ['2', '0', '2', '3', '4', '1'],
             ['2', '0', '2', '4', '1', '3'],
             ['2', '0', '2', '4', '3', '1'],
             ['2', '0', '3', '1', '2', '4'],
             ['2', '0', '3', '1', '4', '2'],
             ['2', '0', '3', '2', '1', '4'],
             ['2', '0', '3', '2', '4', '1'],
             ['2', '0', '3', '4', '1', '2'],
             ['2', '0', '3', '4', '2', '1'],
             ['2', '0', '4', '1', '2', '3'],
             ['2', '0', '4', '1', '3', '2'],
             ['2', '0', '4', '2', '1', '3'],
             ['2', '0', '4', '2', '3', '1'],
             ['2', '0', '4', '3', '1', '2'],
             ['2', '0', '4', '3', '2', '1'],
             ['2', '1', '0', '2', '3', '4'],
             ['2', '1', '0', '2', '4', '3'],
             ['2', '1', '0', '3', '2', '4'],
             ['2', '1', '0', '3', '4', '2'],
             ['2', '1', '0', '4', '2', '3'],
             ['2', '1', '0', '4', '3', '2'],
             ['2', '1', '2', '0', '3', '4'],
             ['2', '1', '2', '0', '4', '3'],
             ['2', '1', '2', '3', '0', '4'],
             ['2', '1', '2', '3', '4', '0'],
             ['2', '1', '2', '4', '0', '3'],
             ['2', '1', '2', '4', '3', '0'],
             ['2', '1', '3', '0', '2', '4'],
             ['2', '1', '3', '0', '4', '2'],
             ['2', '1', '3', '2', '0', '4'],
             ['2', '1', '3', '2', '4', '0'],
             ['2', '1', '3', '4', '0', '2'],
             ['2', '1', '3', '4', '2', '0'],
             ['2', '1', '4', '0', '2', '3'],
             ['2', '1', '4', '0', '3', '2'],
             ['2', '1', '4', '2', '0', '3'],
             ['2', '1', '4', '2', '3', '0'],
             ['2', '1', '4', '3', '0', '2'],
             ['2', '1', '4', '3', '2', '0'],
             ['2', '2', '0', '1', '3', '4'],
             ['2', '2', '0', '1', '4', '3'],
             ['2', '2', '0', '3', '1', '4'],
             ['2', '2', '0', '3', '4', '1'],
             ['2', '2', '0', '4', '1', '3'],
             ['2', '2', '0', '4', '3', '1'],
             ['2', '2', '1', '0', '3', '4'],
             ['2', '2', '1', '0', '4', '3'],
             ['2', '2', '1', '3', '0', '4'],
             ['2', '2', '1', '3', '4', '0'],
             ['2', '2', '1', '4', '0', '3'],
             ['2', '2', '1', '4', '3', '0'],
             ['2', '2', '3', '0', '1', '4'],
             ['2', '2', '3', '0', '4', '1'],
             ['2', '2', '3', '1', '0', '4'],
             ['2', '2', '3', '1', '4', '0'],
             ['2', '2', '3', '4', '0', '1'],
             ['2', '2', '3', '4', '1', '0'],
             ['2', '2', '4', '0', '1', '3'],
             ['2', '2', '4', '0', '3', '1'],
             ['2', '2', '4', '1', '0', '3'],
             ['2', '2', '4', '1', '3', '0'],
             ['2', '2', '4', '3', '0', '1'],
             ['2', '2', '4', '3', '1', '0'],
             ['2', '3', '0', '1', '2', '4'],
             ['2', '3', '0', '1', '4', '2'],
             ['2', '3', '0', '2', '1', '4'],
             ['2', '3', '0', '2', '4', '1'],
             ['2', '3', '0', '4', '1', '2'],
             ['2', '3', '0', '4', '2', '1'],
             ['2', '3', '1', '0', '2', '4'],
             ['2', '3', '1', '0', '4', '2'],
             ['2', '3', '1', '2', '0', '4'],
             ['2', '3', '1', '2', '4', '0'],
             ['2', '3', '1', '4', '0', '2'],
             ['2', '3', '1', '4', '2', '0'],
             ['2', '3', '2', '0', '1', '4'],
             ['2', '3', '2', '0', '4', '1'],
             ['2', '3', '2', '1', '0', '4'],
             ['2', '3', '2', '1', '4', '0'],
             ['2', '3', '2', '4', '0', '1'],
             ['2', '3', '2', '4', '1', '0'],
             ['2', '3', '4', '0', '1', '2'],
             ['2', '3', '4', '0', '2', '1'],
             ['2', '3', '4', '1', '0', '2'],
             ['2', '3', '4', '1', '2', '0'],
             ['2', '3', '4', '2', '0', '1'],
             ['2', '3', '4', '2', '1', '0'],
             ['2', '4', '0', '1', '2', '3'],
             ['2', '4', '0', '1', '3', '2'],
             ['2', '4', '0', '2', '1', '3'],
             ['2', '4', '0', '2', '3', '1'],
             ['2', '4', '0', '3', '1', '2'],
             ['2', '4', '0', '3', '2', '1'],
             ['2', '4', '1', '0', '2', '3'],
             ['2', '4', '1', '0', '3', '2'],
             ['2', '4', '1', '2', '0', '3'],
             ['2', '4', '1', '2', '3', '0'],
             ['2', '4', '1', '3', '0', '2'],
             ['2', '4', '1', '3', '2', '0'],
             ['2', '4', '2', '0', '1', '3'],
             ['2', '4', '2', '0', '3', '1'],
             ['2', '4', '2', '1', '0', '3'],
             ['2', '4', '2', '1', '3', '0'],
             ['2', '4', '2', '3', '0', '1'],
             ['2', '4', '2', '3', '1', '0'],
             ['2', '4', '3', '0', '1', '2'],
             ['2', '4', '3', '0', '2', '1'],
             ['2', '4', '3', '1', '0', '2'],
             ['2', '4', '3', '1', '2', '0'],
             ['2', '4', '3', '2', '0', '1'],
             ['2', '4', '3', '2', '1', '0'],
             ['3', '0', '1', '2', '2', '4'],
             ['3', '0', '1', '2', '4', '2'],
             ['3', '0', '1', '4', '2', '2'],
             ['3', '0', '2', '1', '2', '4'],
             ['3', '0', '2', '1', '4', '2'],
             ['3', '0', '2', '2', '1', '4'],
             ['3', '0', '2', '2', '4', '1'],
             ['3', '0', '2', '4', '1', '2'],
             ['3', '0', '2', '4', '2', '1'],
             ['3', '0', '4', '1', '2', '2'],
             ['3', '0', '4', '2', '1', '2'],
             ['3', '0', '4', '2', '2', '1'],
             ['3', '1', '0', '2', '2', '4'],
             ['3', '1', '0', '2', '4', '2'],
             ['3', '1', '0', '4', '2', '2'],
             ['3', '1', '2', '0', '2', '4'],
             ['3', '1', '2', '0', '4', '2'],
             ['3', '1', '2', '2', '0', '4'],
             ['3', '1', '2', '2', '4', '0'],
             ['3', '1', '2', '4', '0', '2'],
             ['3', '1', '2', '4', '2', '0'],
             ['3', '1', '4', '0', '2', '2'],
             ['3', '1', '4', '2', '0', '2'],
             ['3', '1', '4', '2', '2', '0'],
             ['3', '2', '0', '1', '2', '4'],
             ['3', '2', '0', '1', '4', '2'],
             ['3', '2', '0', '2', '1', '4'],
             ['3', '2', '0', '2', '4', '1'],
             ['3', '2', '0', '4', '1', '2'],
             ['3', '2', '0', '4', '2', '1'],
             ['3', '2', '1', '0', '2', '4'],
             ['3', '2', '1', '0', '4', '2'],
             ['3', '2', '1', '2', '0', '4'],
             ['3', '2', '1', '2', '4', '0'],
             ['3', '2', '1', '4', '0', '2'],
             ['3', '2', '1', '4', '2', '0'],
             ['3', '2', '2', '0', '1', '4'],
             ['3', '2', '2', '0', '4', '1'],
             ['3', '2', '2', '1', '0', '4'],
             ['3', '2', '2', '1', '4', '0'],
             ['3', '2', '2', '4', '0', '1'],
             ['3', '2', '2', '4', '1', '0'],
             ['3', '2', '4', '0', '1', '2'],
             ['3', '2', '4', '0', '2', '1'],
             ['3', '2', '4', '1', '0', '2'],
             ['3', '2', '4', '1', '2', '0'],
             ['3', '2', '4', '2', '0', '1'],
             ['3', '2', '4', '2', '1', '0'],
             ['3', '4', '0', '1', '2', '2'],
             ['3', '4', '0', '2', '1', '2'],
             ['3', '4', '0', '2', '2', '1'],
             ['3', '4', '1', '0', '2', '2'],
             ['3', '4', '1', '2', '0', '2'],
             ['3', '4', '1', '2', '2', '0'],
             ['3', '4', '2', '0', '1', '2'],
             ['3', '4', '2', '0', '2', '1'],
             ['3', '4', '2', '1', '0', '2'],
             ['3', '4', '2', '1', '2', '0'],
             ['3', '4', '2', '2', '0', '1'],
             ['3', '4', '2', '2', '1', '0'],
             ['4', '0', '1', '2', '2', '3'],
             ['4', '0', '1', '2', '3', '2'],
             ['4', '0', '1', '3', '2', '2'],
             ['4', '0', '2', '1', '2', '3'],
             ['4', '0', '2', '1', '3', '2'],
             ['4', '0', '2', '2', '1', '3'],
             ['4', '0', '2', '2', '3', '1'],
             ['4', '0', '2', '3', '1', '2'],
             ['4', '0', '2', '3', '2', '1'],
             ['4', '0', '3', '1', '2', '2'],
             ['4', '0', '3', '2', '1', '2'],
             ['4', '0', '3', '2', '2', '1'],
             ['4', '1', '0', '2', '2', '3'],
             ['4', '1', '0', '2', '3', '2'],
             ['4', '1', '0', '3', '2', '2'],
             ['4', '1', '2', '0', '2', '3'],
             ['4', '1', '2', '0', '3', '2'],
             ['4', '1', '2', '2', '0', '3'],
             ['4', '1', '2', '2', '3', '0'],
             ['4', '1', '2', '3', '0', '2'],
             ['4', '1', '2', '3', '2', '0'],
             ['4', '1', '3', '0', '2', '2'],
             ['4', '1', '3', '2', '0', '2'],
             ['4', '1', '3', '2', '2', '0'],
             ['4', '2', '0', '1', '2', '3'],
             ['4', '2', '0', '1', '3', '2'],
             ['4', '2', '0', '2', '1', '3'],
             ['4', '2', '0', '2', '3', '1'],
             ['4', '2', '0', '3', '1', '2'],
             ['4', '2', '0', '3', '2', '1'],
             ['4', '2', '1', '0', '2', '3'],
             ['4', '2', '1', '0', '3', '2'],
             ['4', '2', '1', '2', '0', '3'],
             ['4', '2', '1', '2', '3', '0'],
             ['4', '2', '1', '3', '0', '2'],
             ['4', '2', '1', '3', '2', '0'],
             ['4', '2', '2', '0', '1', '3'],
             ['4', '2', '2', '0', '3', '1'],
             ['4', '2', '2', '1', '0', '3'],
             ['4', '2', '2', '1', '3', '0'],
             ['4', '2', '2', '3', '0', '1'],
             ['4', '2', '2', '3', '1', '0'],
             ['4', '2', '3', '0', '1', '2'],
             ['4', '2', '3', '0', '2', '1'],
             ['4', '2', '3', '1', '0', '2'],
             ['4', '2', '3', '1', '2', '0'],
             ['4', '2', '3', '2', '0', '1'],
             ['4', '2', '3', '2', '1', '0'],
             ['4', '3', '0', '1', '2', '2'],
             ['4', '3', '0', '2', '1', '2'],
             ['4', '3', '0', '2', '2', '1'],
             ['4', '3', '1', '0', '2', '2'],
             ['4', '3', '1', '2', '0', '2'],
             ['4', '3', '1', '2', '2', '0'],
             ['4', '3', '2', '0', '1', '2'],
             ['4', '3', '2', '0', '2', '1'],
             ['4', '3', '2', '1', '0', '2'],
             ['4', '3', '2', '1', '2', '0'],
             ['4', '3', '2', '2', '0', '1'],
             ['4', '3', '2', '2', '1', '0']]
    for i, j in izip_longest(util.next_permutations("012234"),
                             valid[1:], fillvalue='0'):
        testing.assert_array_equal(i, j)
    for i, j in izip_longest(util.next_permutations(valid[14]),
                             valid[15:], fillvalue='0'):
        testing.assert_array_equal(i, j)
    for i, j in izip_longest(util.next_permutations(iter(("0", "1", "2",
                                                          "2", "3", "4"))),
                             valid[1:], fillvalue='0'):
        testing.assert_array_equal(i, j)
    for i, j in izip_longest(util.next_permutations(("0", "1", "2",
                                                     "2", "3", "4")),
                             valid[1:], fillvalue='0'):
        testing.assert_array_equal(i, j)
    for i in util.next_permutations(32541):
        raise AssertionError("Error: next_permutation failed for int(32541).")
    for i in util.next_permutations("0"):
        raise AssertionError("Error: next_permutation failed for '0'.")
    for i in util.next_permutations("000"):
        raise AssertionError("Error: next_permutation failed for '000'.")


def test_probs_to_class():
    arr = "42"
    testing.assert_raises(TypeError, util.probs_to_class, arr)
    arr = np.asarray([[0.1, 0.5, 0.4], [0.2, 0.8, 0.0], [0.1, 0.45, 0.45]])
    testing.assert_array_equal(util.probs_to_class(arr),
                               np.asarray([1, 1, 1]))
    arr = [np.asarray([[0.1, 0.5, 0.4], [0.2, 0.8, 0.0], [0.1, 0.45, 0.45]])]
    testing.assert_array_equal(util.probs_to_class(arr),
                               np.asarray([1, 1, 1]))
    arr = [np.asarray([[0.1, 0.9], [0.2, 0.8], [0.5, 0.5]]),
           np.asarray([[0.8, 0.2], [0.3, 0.7], [0.0, 1.0]])]
    testing.assert_array_equal(np.asarray([[1, 0], [1, 1], [0, 1]]),
                               util.probs_to_class(arr))


def test_round_nearest():
    test = np.random.rand(30, 20)
    testcopy = np.copy(test)
    nearest = np.copy(test)
    for x in np.nditer(nearest, op_flags=['readwrite']):
        x[...] = 0 if (x < 0.5) else 1
    testing.assert_array_equal(
        util.round_nearest(testcopy, inplace=False, threshold=0.5), nearest)
    try:
        testing.assert_array_equal(testcopy, nearest)
    except AssertionError:
        pass
    else:
        testing.assert_true(
            False, msg="Error: round_nearest not inplace not working.")
    testing.assert_array_equal(
        util.round_nearest(testcopy, inplace=True, threshold=0.5), nearest)
    testing.assert_array_equal(testcopy, nearest)
    testcopy = np.copy(test)
    nearest = np.copy(test)
    for x in np.nditer(nearest, op_flags=['readwrite']):
        x[...] = 0 if (x < 0.3) else 1
    testing.assert_array_equal(
        util.round_nearest(testcopy, inplace=False, threshold=0.3), nearest)


def test_shuffle_array():
    random_state = 0
    original_array = np.array(list(xrange(0, 100)))
    array = original_array.copy()
    check_random_state(random_state).shuffle(array)
    testing.assert_array_equal(util.shuffle_array(
        original_array, inplace=False, random_state=random_state), array)
    testing.assert_array_equal(original_array, np.array(list(xrange(0, 100))))
    testing.assert_array_equal(util.shuffle_array(
        original_array, inplace=True, random_state=random_state), array)
    testing.assert_array_equal(original_array, array)
    util.shuffle_array(original_array, inplace=True,
                       random_state=random_state)
    check_random_state(random_state + 1).shuffle(array)
    testing.assert_equal(len(array), len(original_array))
    testing.assert_true(not np.alltrue(array == original_array))
