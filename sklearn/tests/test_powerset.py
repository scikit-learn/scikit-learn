# -*- coding: utf8 -*-
# Author: Alain Pena <alain.pena.r@gmail.com>
#
# License: BSD 3 clause

import numpy as np

from sklearn.utils import testing
from sklearn.datasets import make_multilabel_classification
from sklearn.tree import DecisionTreeClassifier
from sklearn.powerset import Powerset


def test_identities():
    # init
    random_state = 0
    n_features = 10
    n_labels = 2
    n_samples = 10
    n_classes = 5
    length = 1
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
    y_predict = tree.predict(Xt)
    y_probas = tree.predict_proba(Xt)
    tree.fit(X, y[:, 0].ravel())
    y_probasr = tree.predict_proba(Xt)

    # Checks uninitialized powerset
    power = Powerset()
    testing.assert_raises(ValueError, power.probas_unpowerize, y_predict)
    testing.assert_raises(ValueError, power.unpowerize, y_predict)
    testing.assert_raises(ValueError, power.majority_unpowerize, y_predict)

    # amplified un/powerize methods
    amplified = power.unpowerize(power.amplified_powerize(yt))
    testing.assert_array_equal(amplified, yt)
    amplified = power.unpowerize(power.amplified_powerize(yt[0, :].ravel()))
    testing.assert_array_equal(amplified, yt[0, :].ravel())
    power._Powerset__amplified_unpowerize(amplified)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__compressed_unpowerize, amplified)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__compressed_string_unpowerize, amplified)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__dummy_unpowerize, amplified)
    tree.fit(X, power.amplified_powerize(y))
    amplified = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probas, amplified, atol=1)
    tree.fit(X, power.amplified_powerize(y[:, 0].ravel()))
    amplified = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probasr, amplified, atol=1)
    # compressed un/powerize methods
    compressed = power.unpowerize(power.compressed_powerize(yt))
    testing.assert_array_equal(compressed, yt)
    compressed = power.unpowerize(power.compressed_powerize(yt[0, :].ravel()))
    testing.assert_array_equal(compressed, yt[0, :].ravel())
    compressed = power.unpowerize(power.compressed_powerize(yt[0, :]))
    testing.assert_array_equal(compressed, yt[0, :])
    compressed = power.unpowerize(power.compressed_powerize(yt[:, 0]))
    testing.assert_array_equal(compressed, yt[:, 0])
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__amplified_unpowerize, compressed)
    power._Powerset__compressed_unpowerize(compressed)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__compressed_string_unpowerize, compressed)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__dummy_unpowerize, compressed)
    tree.fit(X, power.compressed_powerize(y))
    compressed = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probas, compressed, atol=1)
    tree.fit(X, power.compressed_powerize(y[:, 0].ravel()))
    compressed = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probasr, compressed, atol=1)  # --
    # dummy un/powerize methods
    dummy = power.unpowerize(power.dummy_powerize(yt))
    testing.assert_array_equal(dummy, yt)
    dummy = power.unpowerize(power.dummy_powerize(yt[0, :].ravel()))
    testing.assert_array_equal(dummy, yt[0, :].ravel())
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__amplified_unpowerize, dummy)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__compressed_unpowerize, dummy)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__compressed_string_unpowerize, dummy)
    power._Powerset__dummy_unpowerize(dummy)
    tree.fit(X, power.dummy_powerize(y))
    dummy = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probas, dummy, atol=1)
    tree.fit(X, power.dummy_powerize(y[:, 0].ravel()))
    dummy = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probasr, dummy, atol=1)
    # compressed string un/powerize methods
    compressed_string = power.unpowerize(power.compressed_string_powerize(yt))
    testing.assert_array_equal(compressed_string, yt)
    testing.assert_array_equal(power.compressed_string_powerize(yt),
                               power.compressed_powerize(yt))
    compressed_string = power.unpowerize(
        power.compressed_string_powerize(yt[0, :].ravel()))
    testing.assert_array_equal(compressed_string, yt[0, :].ravel())
    testing.assert_array_equal(power.compressed_powerize(yt[0, :].ravel()),
                               power.compressed_string_powerize(
                                   yt[0, :].ravel())
                               )
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__amplified_unpowerize, compressed_string)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__compressed_unpowerize, compressed_string)
    power._Powerset__compressed_string_unpowerize(compressed_string)
    testing.assert_raise_message(
        ValueError, 'Wrong unpowerizing method!',
        power._Powerset__dummy_unpowerize, compressed_string)
    tree.fit(X, power.compressed_string_powerize(y))
    compressed_string = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probas, compressed_string, atol=1)
    tree.fit(X, power.compressed_string_powerize(y[:, 0].ravel()))
    compressed_string = power.probas_unpowerize(tree.predict_proba(Xt))
    testing.assert_allclose(y_probasr, compressed_string, atol=1)


def test_doc():
    powerset = Powerset()
    amplified = [[1, 0, 1],
                 [0, 1, 1],
                 [0, 0, 0]]
    amplified = Powerset().amplified_powerize(np.array(amplified))
    testing.assert_array_equal(amplified, np.array([5, 3, 0]))
    compressed_string = [[1, 0, 1],
                         [0, 1, 1],
                         [0, 0, 0]]
    compressed_string = Powerset().compressed_string_powerize(
        np.array(compressed_string))
    testing.assert_array_equal(compressed_string, np.array([0, 1, 2]))
    compressed = [[1, 0, 1],
                  [0, 1, 1],
                  [0, 0, 0]]
    compressed = Powerset().compressed_powerize(np.array(compressed))
    testing.assert_array_equal(compressed, np.array([0, 1, 2]))
    dummy = np.array([[1, 0, 1],
                      [0, 1, 1],
                      [0, 0, 0]])
    dummy_ = Powerset().dummy_powerize(dummy)
    testing.assert_array_equal(dummy, dummy_)
    p = powerset.amplified_powerize(np.array(
        [[0, 0, 0, ], [0, 0, 1, ], [0, 1, 0, ], [0, 1, 1, ],
         [1, 0, 0, ], [1, 0, 1, ], [1, 1, 0, ], [1, 1, 1, ]]))
    testing.assert_array_equal(p, [0, 1, 2, 3, 4, 5, 6, 7])
    c = np.array([[0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16]])
    r = [np.array([[0.42, 0.58]]),
         np.array([[0.46, 0.54]]),
         np.array([[0.48, 0.52]])]
    retv = powerset.probas_unpowerize(c)
    testing.assert_allclose(r, retv, rtol=1e-10)

    p = powerset.amplified_powerize(np.array(
        [[0, 0, 0, ], [0, 0, 1, ], [0, 1, 0, ], [0, 1, 1, ],
         [1, 0, 0, ], [1, 0, 1, ], [1, 1, 0, ], [1, 1, 1, ]]))
    testing.assert_array_equal(p, [0, 1, 2, 3, 4, 5, 6, 7])
    c = np.array([[0.09, 0.10, 0.11, 0.12, 0.13, 0.16, 0.15, 0.14],
                  [0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16]])
    r = [np.array([1, 1]),
         np.array([0, 1]),
         np.array([1, 1])]
    retv = powerset.majority_unpowerize(c)
    testing.assert_allclose(r, retv, rtol=1e-10)
    array = [[1, 0, 1],
             [0, 1, 1],
             [0, 0, 0]]
    p = powerset.unpowerize(powerset.compressed_powerize(np.array(array)))
    testing.assert_array_equal(p, array)
    array = [[1, 0, 1],
             [0, 1, 1],
             [0, 0, 0]]
    p = powerset.universal_powerize(np.array(array))
    testing.assert_array_equal(p, np.array(["101", "011", "000"]))
    p = powerset.universal_unpowerize(np.array(["101", "011", "000"]))
    testing.assert_array_equal(p, np.array(array))


def test_values():
    test = np.random.random_integers(0, 1, (10, 32))
    test[0, 0] = 1
    powerset = Powerset()
    amplified = []

    def _ampl(t):
        amplified[:] = []
        amplified.append(powerset.amplified_powerize(t))

    try:
        long(5)
    except NameError:
        _ampl(test)
    else:
        try:
            testing.assert_warns_message(UserWarning, "Overflow", _ampl, test)
        except AssertionError:
            test = np.random.random_integers(0, 1, (10, 64))
            test[0, 0] = 1
            testing.assert_warns_message(UserWarning, "Overflow", _ampl, test)

    amplified = amplified[0]
    amplified_un = powerset.unpowerize(amplified)

    compressed = powerset.compressed_powerize(test)
    uncompressed = powerset.unpowerize(compressed)

    testing.assert_equal(test.shape, amplified_un.shape)
    testing.assert_array_equal(test, amplified_un)
    testing.assert_equal(test.shape, uncompressed.shape)
    testing.assert_array_equal(test, uncompressed)
    testing.assert_equal(compressed.shape, amplified.shape)
    testing.assert_equal(test.shape[0], compressed.shape[0])

    # Mapping: class 0 means the labels [0,1,0,0]
    u = {0: [0, 1, 0, 0], 1: [0, 0, 1, 1], 2: [1, 1, 0, 0]}
    # Result of predict_proba
    c = [0.4, 0.5, 0.1]
    # Corrected result
    rg = np.array([0, 0, 1, 0])
    r = np.array([[0.4 + 0.5, 0.1],
                 [0.5, 0.4 + 0.1],
                 [0.4 + 0.1, 0.5],
                 [0.4 + 0.1, 0.5]])
    r2 = [[r[0], r[0]], [r[1], r[1]], [r[2], r[2]], [r[3], r[3]]]
    r3 = [np.array(x) for x in r2]
    retv = Powerset._probas_aux(uncompressed=u, cls_probas_array=c)
    testing.assert_array_equal(r, retv)
    ur = np.array([[0, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 1], [1, 1, 0, 0]])
    p = powerset.compressed_powerize(ur)
    testing.assert_array_equal(p, [0, 0, 1, 2])
    cm = np.empty((2, len(c)), dtype=float)
    cm[0, :] = np.asarray(c)
    cm[1, :] = np.asarray(c)
    cm_ = [cm]
    cms = [[0, 0], [0, 0], [1, 1], [1, 1]]
    cm2 = np.empty((4, len(c)), dtype=float)
    cm2[0, :] = np.asarray(c)
    cm2[1, :] = np.asarray(c)
    cm2[2, :] = np.asarray(c)
    cm2[2, 1] = 0.3
    cm2[2, 2] = 1 - cm2[2, 0] - 0.3
    cm2[3, :] = [0.33, 0.33, 0.33]
    cm2s = [[0, 0, 0, 0], [0, 0, 1, 1], [1, 1, 0, 0], [1, 1, 0, 0]]
    p = powerset.probas_unpowerize(cm)
    testing.assert_allclose(p, r2, rtol=1e-10)
    testing.assert_allclose(p, r3, rtol=1e-10)
    p = powerset.majority_unpowerize(cm)
    testing.assert_allclose(p, cms, rtol=1e-10)
    p = powerset.probas_unpowerize(cm_)
    testing.assert_allclose(p, r2, rtol=1e-10)
    testing.assert_allclose(p, r3, rtol=1e-10)
    p = powerset.majority_unpowerize(cm_)
    testing.assert_allclose(p, cms, rtol=1e-10)
    testing.assert_raise_message(
        ValueError, 'single', powerset.majority_unpowerize, [cm, cm])
    p = powerset.majority_unpowerize(cm2)
    testing.assert_allclose(p, cm2s, rtol=1e-10)
    p = powerset.dummy_powerize(r)
    p = powerset.probas_unpowerize(r)
    testing.assert_array_equal(p, r)
    p = powerset.dummy_powerize(rg)
    p = powerset.majority_unpowerize(r)
    testing.assert_array_equal(p, np.argmax(r, axis=1))
    test += 1
    testing.assert_raise_message(
        ValueError, 'binary', powerset.compressed_powerize, test)
    test -= 1
    test[0, 0] = 2
    testing.assert_raise_message(
        ValueError, 'not contain 2 classes',
        powerset.compressed_powerize, test)
    test[0, 0] = 1
    testing.assert_raise_message(
        ValueError, 'bidimensional', powerset.compressed_powerize, r2)
    p = powerset.compressed_string_powerize(ur)
    p = powerset.probas_unpowerize(cm)
    testing.assert_raise_message(ValueError, 'invalid', powerset.unpowerize, p)
    testing.assert_raise_message(
        ValueError, 'dimensional', Powerset.universal_unpowerize, p)
    p = powerset.compressed_string_powerize(ur)
    p = powerset.probas_unpowerize(cm)
    testing.assert_array_equal(p, powerset.probas_unpowerize([cm]))
