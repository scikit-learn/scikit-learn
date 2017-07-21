import numpy as np

from sklearn.glvq.glvq import GlvqModel
from sklearn.glvq.gmlvq import GmlvqModel
from sklearn.glvq.grlvq import GrlvqModel
from sklearn.glvq.lgmlvq import LgmlvqModel
from sklearn.utils.testing import assert_greater, assert_raise_message, assert_equal

from sklearn import datasets
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
    model = GlvqModel() #, display=True)
    model.fit(iris.data, iris.target)
    assert_greater(model.score(iris.data, iris.target), 0.89)
    print(model.score(iris.data, iris.target))

    model = GlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                      initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', GlvqModel(max_iter='5').fit, iris.data, iris.target)

    # TODO: find test


def test_grlvq_iris():
    print('test grlvq')
    model = GrlvqModel(initial_prototypes=[[0,2],[1,6]],initial_rototype_labels=[1,2]) #, display=True)
    #model.fit(iris.data, iris.target)
    X = np.array([[0,0],[0,4],[1,4],[1,8]])
    y = np.array([1,1,2,2])
    model.fit(X,y)
    #assert_greater(model.score(iris.data, iris.target), 0.95)

    model = GrlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                      initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', GrlvqModel(max_iter='5').fit, iris.data, iris.target)


def test_gmlvq_iris():
    print('test gmlvq')
    model = GmlvqModel(initial_prototypes=[[0,2],[1,6]],initial_rototype_labels=[1,2]) #, display=True)
    #model.fit(iris.data, iris.target)
    X = np.array([[0,0],[0,4],[1,4],[1,8]])
    y = np.array([1,1,2,2])
    model.fit(X,y)
    #assert_greater(model.score(iris.data, iris.target), 0.95)
    print(model.score(X,y))

    model = GmlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                      initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', GmlvqModel(max_iter='5').fit, iris.data, iris.target)


def test_lgmlvq_iris():
    print('test lgmlvq')
    model = LgmlvqModel(initial_prototypes=[[0,2],[1,6]],initial_rototype_labels=[1,2],initial_matrices=[np.ones([2,2]),np.ones([2,2])],dim=[2,2]) #, display=True)
    #model.fit(iris.data, iris.target)
    X = np.array([[0,0],[0,4],[1,4],[1,8]])
    y = np.array([1,1,2,2])
    model.fit(X,y)
    #assert_greater(model.score(iris.data, iris.target), 0.9)
    print(model.score(iris.data, iris.target))

    model = LgmlvqModel(gtol=1e-4, max_iter=5, prototypes_per_class=2,
                      initial_prototypes=[[1, 1, 1, 1], [1, 1, 1, 1], [2, 2, 2, 2]], initial_rototype_labels=[0, 1, 2])
    assert_raise_message(ValueError, 'The initial prototypes have wrong shape', model.fit, iris.data, iris.target)
    assert_raise_message(ValueError, 'max_iter must be an integer', LgmlvqModel(max_iter='5').fit, iris.data, iris.target)

    # TODO: find test[[ 0.05592003  0.01807518]

def test():
    g = LgmlvqModel()
    variables = np.array([[-0.0422,0.1131],
                          [3.1333,6.9411],
                          [0.6294,-0.7460],
                          [0.8116,0.8268],
                          [0.2647,-0.4430],
                          [-0.8049,0.0938]])
    training = np.array([[0,0],[0,4],[1,4],[1,8]])
    labels = np.array([[True,False],[True,False],[False,True],[False,True]])
    g.c_w_ = np.array([1,2])
    g.dim_ = [2,2]
    g.regularization = np.repeat(0, 2)
    g.g(variables,training,labels,random_state=np.random.RandomState(),lr_relevances=1,lr_prototypes=0)