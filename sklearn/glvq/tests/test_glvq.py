import numpy as np

from sklearn.glvq.glvq import GlvqModel
from sklearn.glvq.gmlvq import GmlvqModel
from sklearn.glvq.grlvq import GrlvqModel
from sklearn.glvq.lgmlvq import LgmlvqModel
from sklearn.utils.testing import assert_greater, assert_raise_message, assert_allclose

from sklearn import datasets
from sklearn.utils import check_random_state

# also load the iris dataset
iris = datasets.load_iris()
rng = check_random_state(42)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


# TODO: grlvq not working with lbfgs-b

def test_glvq_iris():
    print('test glvq')
    model = GlvqModel()  # , display=True)
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.94)

    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.94)
    print(model.score(iris.data, iris.target))

    assert_raise_message(ValueError, 'display must be a boolean',
                         GlvqModel(display='true').fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'gtol must be a positive float',
                         GlvqModel(gtol=-1.0).fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape',
                         GlvqModel(initial_prototypes=[[1, 1], [2, 2]]).fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'Prototype labels and test data classes dont match',
                         GlvqModel(initial_prototypes=[[1, 1, 1, 1, 'a'], [2, 2, 2, 2, 5], [2, 2, 2, 2, -3]]).fit,
                         iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an positive integer', GlvqModel(max_iter='5').fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'max_iter must be an positive integer', GlvqModel(max_iter=0).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'max_iter must be an positive integer', GlvqModel(max_iter=-1).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'Values in prototypes_per_class must be greater than 0',
                         GlvqModel(prototypes_per_class=np.zeros(np.unique(iris.target).size) - 1).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'Length of prototypes per class does not fit the number of',
                         GlvqModel(prototypes_per_class=[1, 2]).fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'X has wrong number of features', model.predict, [[1, 2], [3, 4]])


def test_grlvq_iris():
    print('test grlvq')
    model = GrlvqModel(initial_prototypes=[[0, 2, 1], [1, 6, 2]])  # , display=True)
    X = np.array([[0, 0], [0, 4], [1, 4], [1, 8]])
    y = np.array([1, 1, 2, 2])
    model.fit(X, y)
    assert_allclose(np.array([1.0, 0.0]), model.lambda_, rtol=0.1, atol=0.1)

    model = GrlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                       initial_prototypes=[[1, 1, 1, 1, 0], [1, 1, 1, 1, 1], [2, 2, 2, 2, 2]])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an positive integer', GrlvqModel(max_iter='5').fit, iris.data,
                         iris.target)


def test_gmlvq_iris():
    print('test gmlvq')
    model = GmlvqModel(initial_prototypes=[[0, 2], [1, 6]], initial_rototype_labels=[1, 2])  # , display=True)
    X = np.array([[0, 0], [0, 4], [1, 4], [1, 8]])
    y = np.array([1, 1, 2, 2])
    model.fit(X, y)
    assert_allclose(np.array([[0.7, 0.4], [0.4, 0.1]]), model.omega_, rtol=0.2)
    print(model.omega_)
    print(model.score(X, y))

    model = GmlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                       initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', GmlvqModel(max_iter='5').fit, iris.data,
                         iris.target)


def test_lgmlvq_iris():
    print('test lgmlvq')
    model = LgmlvqModel(initial_prototypes=[[0, 2], [1, 6]], initial_rototype_labels=[1, 2],
                        initial_matrices=[np.ones([2, 2]), np.ones([2, 2])], dim=[2, 2], display=True)
    X = np.array([[0, 0], [0, 4], [1, 4], [1, 8]])
    y = np.array([1, 1, 2, 2])
    model.fit(X, y)
    print(model.score(X, y))

    model = LgmlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                        initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]],
                        initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', LgmlvqModel(max_iter='5').fit, iris.data,
                         iris.target)
