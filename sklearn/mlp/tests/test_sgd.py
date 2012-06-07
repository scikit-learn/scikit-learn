import numpy as np

from ..classes import MLPClassifier
from ..classes_andreas import MLPClassifierA
from ... import datasets, preprocessing
from nose.tools import assert_true
from datetime import datetime


def test_cross_entropy():
    X, y = datasets.make_classification()

    X = preprocessing.scale(X)

    classifier = MLPClassifier(n_hidden=10, lr=0.3, batch_size=100)
    classifier.fit(X, y,
        max_epochs=500)

    assert_true(np.mean(classifier.predict(X) == y) > .95)


def test_numbers():
    digits = datasets.load_digits()
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))

    data = preprocessing.scale(data)

    classifier = MLPClassifier(n_hidden=10, lr=0.3, batch_size=100)
    classifier.fit(data, digits.target,
        max_epochs=500)

    assert_true(np.mean(classifier.predict(data) == digits.target) > .95)


def test_iris():
    iris = datasets.load_iris()
    data, target = iris.data, iris.target
    data = preprocessing.scale(data)

    classifier = MLPClassifier(n_hidden=10, lr=0.3, batch_size=1)
    classifier.fit(data, target, max_epochs=300)

    assert_true(np.mean(classifier.predict(data) == target) > .95)


def test_sgd_faster():
    digits = datasets.load_digits()
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))

    data = preprocessing.scale(data)
    classifier = MLPClassifierA(n_hidden=10, lr=0.3, loss='cross_entropy', output_layer='softmax', batch_size=100)
    start = datetime.now()
    np.random.seed(0)
    classifier.fit(data[:n_samples / 2], digits.target[:n_samples / 2],
        max_epochs=1000, shuffle_data=False)
    time_sgd_a = datetime.now() - start
    expected_a = digits.target[n_samples / 2:]
    predicted_a = classifier.predict(data[n_samples / 2:])

    classifier = MLPClassifier(n_hidden=10, lr=0.3, loss_function='cross-entropy', output_function='softmax', batch_size=100)
    start = datetime.now()
    np.random.seed(0)
    classifier.fit(data[:n_samples / 2], digits.target[:n_samples / 2],
        max_epochs=1000)
    time_sgd = datetime.now() - start
    expected = digits.target[n_samples / 2:]
    predicted = classifier.predict(data[n_samples / 2:])

    print "SGD algorithm: %s (%f)" % (time_sgd, np.mean(expected == predicted))
    print "Pure implementation: %s (%f)" % (time_sgd_a, np.mean(expected_a == predicted_a))
