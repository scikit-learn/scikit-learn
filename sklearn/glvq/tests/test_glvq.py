import numpy as np

from sklearn.glvq.glvq import GlvqModel
from sklearn.glvq.grlvq import GrlvqModel
from sklearn.glvq.gmlvq import GmlvqModel
from sklearn.glvq.lgmlvq import LgmlvqModel
from sklearn.utils.testing import assert_greater, assert_raise_message, \
    assert_allclose

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
    model = GlvqModel()
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.94)

    assert_raise_message(ValueError, 'display must be a boolean',
                         GlvqModel(display='true').fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'gtol must be a positive float',
                         GlvqModel(gtol=-1.0).fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'the initial prototypes have wrong shape',
                         GlvqModel(initial_prototypes=[[1, 1], [2, 2]]).fit,
                         iris.data, iris.target)
    assert_raise_message(ValueError,
                         'prototype labels and test data classes do not match',
                         GlvqModel(initial_prototypes=[[1, 1, 1, 1, 'a'],
                                                       [2, 2, 2, 2, 5],
                                                       [2, 2, 2, 2, -3]]).fit,
                         iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an positive integer',
                         GlvqModel(max_iter='5').fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'max_iter must be an positive integer',
                         GlvqModel(max_iter=0).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'max_iter must be an positive integer',
                         GlvqModel(max_iter=-1).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError,
                         'values in prototypes_per_class must be positive',
                         GlvqModel(prototypes_per_class=np.zeros(
                             np.unique(iris.target).size) - 1).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError,
                         'length of prototypes per class'
                         ' does not fit the number of',
                         GlvqModel(prototypes_per_class=[1, 2]).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'X has wrong number of features',
                         model.predict, [[1, 2], [3, 4]])


def test_grlvq_iris():
    model = GrlvqModel(regularization=0.5)
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.89)

    model = GrlvqModel(initial_prototypes=[[0, 0, 0], [4, 4, 1]])
    nb_ppc = 10
    X = np.append(
        np.random.multivariate_normal([0, 0], np.array([[0.3, 0], [0, 4]]),
                                      size=nb_ppc),
        np.random.multivariate_normal([4, 4], np.array([[0.3, 0], [0, 4]]),
                                      size=nb_ppc), axis=0)
    y = np.append(np.zeros(nb_ppc), np.ones(nb_ppc), axis=0)
    model.fit(X, y)
    assert_allclose(np.array([1.0, 0.0]), model.lambda_, atol=0.2)

    assert_raise_message(ValueError, 'length of initial relevances is wrong',
                         GrlvqModel(initial_relevances=[1, 2]).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'regularization must be a positive float',
                         GrlvqModel(regularization=-1.0).fit, iris.data,
                         iris.target)


def test_gmlvq_iris():
    model = GmlvqModel(regularization=0.5)
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.94)

    model = GmlvqModel(initial_prototypes=[[0, 0, 0], [4, 4, 1]])
    nb_ppc = 10
    X = np.append(
        np.random.multivariate_normal([0, 0], np.array([[0.3, 0], [0, 4]]),
                                      size=nb_ppc),
        np.random.multivariate_normal([4, 4], np.array([[0.3, 0], [0, 4]]),
                                      size=nb_ppc), axis=0)
    y = np.append(np.zeros(nb_ppc), np.ones(nb_ppc), axis=0)
    model.fit(X, y)
    assert_allclose(np.array([[1, 0], [0.2, 0]]), model.omega_, atol=0.3)

    assert_raise_message(ValueError, 'regularization must be a positive float',
                         GmlvqModel(regularization=-1.0).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError,
                         'initial matrix has wrong number of features',
                         GmlvqModel(
                             initial_matrix=[[1, 2], [3, 4], [5, 6]]).fit,
                         iris.data, iris.target)
    assert_raise_message(ValueError, 'dim must be an positive int',
                         GmlvqModel(dim=0).fit, iris.data, iris.target)


def test_lgmlvq_iris():
    model = LgmlvqModel()
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.95)

    assert_raise_message(ValueError, 'regularization must be a positive float',
                         LgmlvqModel(regularization=-1.0).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError,
                         'length of regularization'
                         ' must be number of prototypes',
                         LgmlvqModel(regularization=[-1.0]).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError,
                         'length of regularization must be number of classes',
                         LgmlvqModel(regularization=[-1.0],
                                     classwise=True).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'initial matrices must be a list',
                         LgmlvqModel(initial_matrices=np.array(
                             [[1, 2], [3, 4], [5, 6]])).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'length of matrices wrong',
                         LgmlvqModel(
                             initial_matrices=[[[1, 2], [3, 4], [5, 6]]]).fit,
                         iris.data, iris.target)
    assert_raise_message(ValueError, 'each matrix should have',
                         LgmlvqModel(
                             initial_matrices=[[[1]], [[1]], [[1]]]).fit,
                         iris.data, iris.target)
    assert_raise_message(ValueError, 'length of matrices wrong',
                         LgmlvqModel(initial_matrices=[[[1, 2, 3]]],
                                     classwise=True).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'each matrix should have',
                         LgmlvqModel(initial_matrices=[[[1]], [[1]], [[1]]],
                                     classwise=True).fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'classwise must be a boolean',
                         LgmlvqModel(classwise="a").fit, iris.data,
                         iris.target)
    assert_raise_message(ValueError, 'dim must be a list of positive ints',
                         LgmlvqModel(dim=[-1]).fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'dim length must be number of prototypes',
                         LgmlvqModel(dim=[1, 1]).fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'dim length must be number of classes',
                         LgmlvqModel(dim=[1, 1], classwise=True).fit,
                         iris.data, iris.target)
    LgmlvqModel(classwise=True, dim=[1], prototypes_per_class=2).fit(
        iris.data, iris.target)

    model = LgmlvqModel(regularization=0.1)
    model.fit(iris.data, iris.target)

    model = LgmlvqModel(initial_prototypes=[[0, 2, 1], [1, 6, 2]],
                        initial_matrices=[np.ones([2, 2]), np.ones([2, 2])],
                        dim=[2, 2])
    X = np.array([[0, 0], [0, 4], [1, 4], [1, 8]])
    y = np.array([1, 1, 2, 2])
    model.fit(X, y)
