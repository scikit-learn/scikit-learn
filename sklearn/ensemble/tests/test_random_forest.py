"""
Testing for Random Forest module (sklearn.ensemble)

"""

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal, \
                          assert_almost_equal
from nose.tools import assert_raises

from sklearn import ensemble, datasets, metrics

# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
Y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]

# also load the iris dataset
# and randomly permute it
iris = datasets.load_iris()
np.random.seed([1]) 
perm = np.random.permutation(iris.target.size)
iris.data = iris.data[perm]
iris.target = iris.target[perm]

# also load the boston dataset
# and randomly permute it
boston = datasets.load_boston()
perm = np.random.permutation(boston.target.size)
boston.data = boston.data[perm]
boston.target = boston.target[perm]

def test_classification_toy():
    """
    Check classification on a toy dataset
    """

    clf = ensemble.RandomForestClassifier(r=0.9)
    clf.fit(X,Y)
    
    assert_array_equal(clf.predict(T), true_result)

    """
    With subsampling
    """
    clf = ensemble.RandomForestClassifier(r=0.9, F=1)
    clf.fit(X,Y)
    
    assert_array_equal(clf.predict(T), true_result)

def test_regression_toy():
    """
    Check regression on a toy dataset
    """
    clf = ensemble.RandomForestRegressor(r=0.9)
    clf.fit(X,Y)

    assert_almost_equal(clf.predict(T), true_result)

    """
    With subsampling
    """
    clf = ensemble.RandomForestRegressor(r=0.9, F=1)
    clf.fit(X,Y)
    
    assert_almost_equal(clf.predict(T), true_result)


def test_iris():
    """
    Check consistency on dataset iris.
    """

    for c in ('gini', \
              'entropy'):
        clf = ensemble.RandomForestClassifier(criterion=c)\
              .fit(iris.data, iris.target)
            
        score = np.mean(clf.predict(iris.data) == iris.target) 
        assert score > 0.9, "Failed with criterion " + c + \
            " and score = " + str(score) 

        clf = ensemble.RandomForestClassifier(criterion=c, F=2)\
              .fit(iris.data, iris.target)
            
        score = np.mean(clf.predict(iris.data) == iris.target) 
        assert score > 0.5, "Failed with criterion " + c + \
            " and score = " + str(score)        

def test_boston():
    """
    Check consistency on dataset boston house prices.
    """
    for c in ('mse',):
        # pass in a seed to ensure that tests using random permutations 
        # are repeatable
        clf = ensemble.RandomForestRegressor(criterion=c, random_state=1)\
              .fit(boston.data, boston.target)
            
        score = np.mean(np.power(clf.predict(boston.data)-boston.target,2)) 
        assert score < 2, "Failed with F=None, criterion " + c + \
            " and score = " + str(score)             

        clf = ensemble.RandomForestRegressor(criterion=c, F=6, random_state=1)\
              .fit(boston.data, boston.target)
        
        #using fewer dimensions reduces the learning ability of this tree, 
        # but reduces training time.
        score = np.mean(np.power(clf.predict(boston.data)-boston.target,2)) 
        assert score < 3, "Failed with F=6, criterion " + c + \
            " and score = " + str(score)    


def test_sanity_checks_predict():
    Xt = np.array(X).T

    clf = ensemble.RandomForestClassifier()
    clf.fit(np.dot(X, Xt), Y)
    assert_raises(ValueError, clf.predict, X)

    clf = ensemble.RandomForestClassifier()
    clf.fit(X, Y)
    assert_raises(ValueError, clf.predict, Xt)



def test_probability():
    """
    Predict probabilities using RandomForestClassifier
    """

    clf = ensemble.RandomForestClassifier()
    clf.fit(iris.data, iris.target)

    prob_predict = clf.predict_proba(iris.data)
    assert_array_almost_equal(
        np.sum(prob_predict, 1), np.ones(iris.data.shape[0]))
    assert np.mean(np.argmax(prob_predict, 1) 
                   == clf.predict(iris.data)) > 0.9

    assert_almost_equal(clf.predict_proba(iris.data),
                        np.exp(clf.predict_log_proba(iris.data)), 8)

def test_error():
    """
    Test that it gives proper exception on deficient input
    """

    #invalid number of trees

if __name__ == '__main__':
    import nose
    nose.runmodule()