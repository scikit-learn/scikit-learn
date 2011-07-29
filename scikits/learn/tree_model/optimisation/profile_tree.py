"""
Profiling for Tree module (scikits.learn.tree_model)

"""

import numpy as np

from scikits.learn import tree_model, datasets, metrics
from scikits.learn.datasets.samples_generator import test_dataset_classif

import cProfile
import pstats

def dtree_classification():

    n_samples = 1000
    dim = 50
    K = 10
    X = np.random.randn(n_samples, dim)
    Y = np.random.randint(0, K, (n_samples,))

    clf = tree_model.DecisionTreeClassifier(K)
    clf.fit(X,Y).predict(X)
    

def dtree_regression():

    n_samples = 1000
    dim = 50
    X = np.random.randn(n_samples, dim)
    Y = np.random.randn(n_samples)

    clf = tree_model.DecisionTreeRegressor()
    clf.fit(X,Y).predict(X)

if __name__ == '__main__':

    cProfile.run('dtree_classification()', 'dtree_classification')
    p = pstats.Stats('dtree_classification')    
    p.sort_stats('time').print_stats()
   
    cProfile.run('dtree_regression()', 'dtree_regression')
    p = pstats.Stats('dtree_regression')    
    p.sort_stats('time').print_stats()