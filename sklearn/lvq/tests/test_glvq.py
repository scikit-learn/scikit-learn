import numpy as np

from sklearn.utils.testing import assert_greater, assert_raise_message

from sklearn import datasets
from sklearn.lvq.glvq import GlvqModel, GrlvqModel, GmlvqModel
from sklearn.utils import check_random_state

# also load the iris dataset
iris = datasets.load_iris()
rng = check_random_state(42)
perm = rng.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]


# TODO: implement grlvq,gmlvq,lgmlvq,lgmrvq

def test_glvq_iris():
    print('test glvq')
    model = GlvqModel(prototypes_per_class=2) #, display=True)
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.89)
    print(model.score(iris.data, iris.target))

    model = GlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                      initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', GlvqModel, max_iter='5')

    # TODO: find test

def test_grlvq_iris():
    print('test grlvq')
    model = GrlvqModel(prototypes_per_class=2) #, display=True)
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.95)
    print(model.score(iris.data, iris.target))

    model = GrlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                      initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', GlvqModel, max_iter='5')

    # TODO: find test

def test_gmlvq_iris():
    print('test gmlvq')
    model = GmlvqModel(prototypes_per_class=2) #, display=True)
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.95)
    print(model.score(iris.data, iris.target))

    model = GmlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                      initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', GlvqModel, max_iter='5')

    # TODO: find test