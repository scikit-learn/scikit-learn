# -*- coding: utf-8 -*-
# Copyright (C) 2008-2011, Luis Pedro Coelho <luis@luispedro.org>
# vim: set ts=4 sts=4 sw=4 expandtab smartindent:
#
# License: MIT. See COPYING.MIT file in the milk distribution
"""
================
Tree Classifier
================

    A decision tree classifier 
    
    Implements Classification and Regression Trees (Breiman et al. 1984)

"""

from __future__ import division
import numpy as np
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ._tree import bincount_k
from ._tree import eval_gini, eval_entropy, eval_miss, eval_mse
import random

__all__ = [
    'DecisionTreeClassifier',
    'DecisionTreeRegressor',
    ]

lookup_c = \
      {'gini' : eval_gini,
       'entropy' : eval_entropy, 
       'miss' : eval_miss,
       }
lookup_r = \
      {'mse' : eval_mse,       
      }   

class Leaf(object):
    '''
        v : target value
            Classification: array-like, shape = [n_features]
                Histogram of target values 
            Regression:  real number
                Mean for the region
    '''
    def __init__(self, v):
        self.v = v
    def __str__(self):
        return 'Leaf(%s)' % (self.v)

class Node(object):
    def __init__(self, featid, featval, left, right):
        self.featid = featid
        self.featval = featval
        self.left = left
        self.right = right
    def __str__(self):
        return 'Node(%s,%s)' % (self.featid, self.featval)
"""    
@TODO Profiling shows that this function is the bottleneck
      now that the intensive evaluation functions have been optimised.
      Consider moving this to _tree.pyx 
"""
def _split(features, labels, criterion):

    n_samples, n_features = features.shape
        
    best = None
    split_error = np.inf
    for i in xrange(n_features):
        domain_i = sorted(set(features[:,i]))
        for d in domain_i[1:]:
            cur_split = (features[:,i] < d)
            error = criterion(labels[cur_split], labels[~cur_split])
            if error < split_error:
                split_error = error
                best = i,d
    return best


def _build_tree(is_classification, features, labels, criterion, \
               max_depth, min_split, F, K):
    
    if len(labels) != len(features):
        raise ValueError("Number of labels does not match number of features\n" +
                         "num labels is %s and num features is %s " % 
                         (len(labels), len(features)))
        
    sample_dims = np.array(xrange(features.shape[1]))    
    if F is not None:
        if F <= 0 :
            raise ValueError("F must be > 0.\n" +
                             "Did you mean to use None to signal no F?")
        if F > features.shape[1] :
            raise ValueError("F must be < num dimensions of features.\n" +
                             "F is %s, n_dims = %s " % (F, features.shape[1]))            
        
        sample_dims = np.sort(np.array(random.sample(xrange(features.shape[1]), F)))
        features = features[:, sample_dims]

    if min_split <= 0:
        raise ValueError("min_split must be greater than zero.\n" + 
                         "min_split is %s." % min_split )
    if max_depth <= 0:
        raise ValueError("max_depth must be greater than zero.\n" + 
                         "max_depth is %s." % max_depth )


    def recursive_partition(features, labels, depth):
        N = float(len(labels))
        if N < min_split or depth >= max_depth:
            if is_classification:
                return Leaf(bincount_k(labels, K) / N ) 
            else:
                return Leaf(np.mean(labels))
        S = _split(features, labels, criterion)
        if S is None:
            if is_classification:
                return Leaf(bincount_k(labels, K) / N ) 
            else:
                return Leaf(np.mean(labels))        
        dim,thresh = S
        split = features[:,dim] < thresh
        return Node(featid=sample_dims[dim],
                    featval=thresh,
                    left =recursive_partition(features[ split], \
                                              labels[ split], depth+1),
                    right=recursive_partition(features[~split], \
                                              labels[~split], depth+1))
        
    return recursive_partition(features, labels, 0)

def _apply_tree(tree, features):
    '''
    conf = apply_tree(tree, features)

    Applies the decision tree to a set of features.
    '''
    if type(tree) is Leaf:
        return tree.v
    if features[tree.featid] < tree.featval:
        return _apply_tree(tree.left, features)
    return _apply_tree(tree.right, features)

"""    
@TODO Print the tree to a .dot file using graphViz
"""
def print_tree(tree):
    '''Print decision tree
    '''
    pass

class BaseDecisionTree(BaseEstimator):
    '''
    Should not be used directly, use derived classes instead
    '''

    _dtree_types = ['classification', 'regression']    
    
    def __init__(self, K, impl, criterion, max_depth, min_split, F, seed):
        
        if not impl in self._dtree_types:
            raise ValueError("impl should be one of %s, %s was given" % (
                self._dtree_types, impl))                    
            
        self.type = impl
        
        if impl == 'classification' and K is None:
            raise ValueError("For classification trees, K (number of classes)\n" +
                             "must be given.")                     
        self.K = K
                
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_split = min_split        
        self.F = F
        
        if seed is not None:
            random.seed(seed)
        
        self.n_features = None
        self.tree = None

    def __str__(self):
        print_tree(self.tree)

    def fit(self, X, y):
        """
        Fit the tree model according to the given training data and
        parameters.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]
            Training vectors, where n_samples is the number of samples
            and n_features is the number of features.

        y : array-like, shape = [n_samples]
            Target values (integers in classification, real numbers in
            regression)
            For classification, labels must correspond to classes 0,1,...,K-1

        Returns
        -------
        self : object
            Returns self.
        """               
        X = np.asanyarray(X, dtype=np.float64, order='C')       
        _, self.n_features = X.shape    
        
        if self.type == 'classification':
            y = np.asanyarray(y, dtype=np.int, order='C') 
            if y.max() >= self.K or y.min() < 0:
                raise ValueError("Labels must be in the range [0 to %s)", 
                                 self.K)  
            self.tree = _build_tree(True, X, y, lookup_c[self.criterion], \
                                    self.max_depth, self.min_split, self.F, \
                                    self.K)
        else: #regression
            y = np.asanyarray(y, dtype=np.float64, order='C')               
            self.tree = _build_tree(False, X, y, lookup_r[self.criterion], \
                                    self.max_depth, self.min_split, self.F, \
                                    None)                  
        return self  
        
    def predict(self, X):
        """
        This function does classification or regression on an array of
        test vectors X.

        For a classification model, the predicted class for each
        sample in X is returned.  For a regression model, the function
        value of X calculated is returned.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples]
        """

        X = np.atleast_2d(X)
        n_samples, n_features = X.shape    
            
        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')
        
        if self.n_features != n_features:
            raise ValueError("Number of features of the model must match the input.\n" + 
                             "Model n_features is %s and input n_features is %s " % \
                             (self.n_features, n_features) )        
        
        C = np.zeros(n_samples, dtype=int)
        for idx, sample in enumerate(X):
            if self.type == 'classification':
                C[idx] = np.argmax(_apply_tree(self.tree, sample))
            else:
                C[idx] = _apply_tree(self.tree, sample)            
        
        return C 


class DecisionTreeClassifier(BaseDecisionTree, ClassifierMixin):
    """Classify a multi-labeled dataset with a decision tree.

    Parameters
    ----------
    K : integer, mandatory 
        number of classes    

    criterion : string
        function to measure goodness of split

    max_depth : integer
        maximum depth of the tree  
              
    min_split : integer
        minimum size to split on
        
    F : integer, optional
        if given, then, choose F features

    seed : integer or array_like, optional
        seed the random number generator
        
    Returns
    -------
    tree : Tree

    
    Example
    -------
    >>> import numpy as np
    >>> from scikits.learn.datasets import load_iris
    >>> from scikits.learn.cross_val import StratifiedKFold
    >>> from scikits.learn import tree_model
    >>> data = load_iris()
    >>> skf = StratifiedKFold(data.target, 10)
    >>> for train_index, test_index in skf:
    ...     tree = tree_model.DecisionTreeClassifier(K=3)
    ...     tree = tree.fit(data.data[train_index], data.target[train_index])
    ...     print np.mean(tree.predict(data.data[test_index]) == data.target[test_index])
    ... 
    0.933333333333
    0.866666666667
    0.8
    0.933333333333
    0.933333333333
    0.933333333333
    0.933333333333
    1.0
    0.866666666667
    1.0

    
    """
    def __init__(self, K, criterion='gini', max_depth=10,\
                  min_split=5, F=None, seed=None):
        BaseDecisionTree.__init__(self, K, 'classification', criterion, \
                                  max_depth, min_split, F, seed)
    
    def predict_proba(self, X):
        """
        This function does classification on a test vector X
        given a model with probability information.

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        P : array-like, shape = [n_samples, n_classes]
            Returns the probability of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.

        """
        X = np.atleast_2d(X)
        n_samples, n_features = X.shape    
            
        if self.tree is None:
            raise Exception('Tree not initialized. Perform a fit first')
        
        if self.n_features != n_features:
            raise ValueError("Number of features of the model must match the input.\n" + 
                             "Model n_features is %s and input n_features is %s " % \
                             (self.n_features, n_features) )        
        
        P = np.zeros((n_samples, self.K))
        for idx, sample in enumerate(X):
            P[idx,:] = _apply_tree(self.tree, sample)
        return P        

    def predict_log_proba(self, X):
        """
        This function does classification on a test vector X

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        P : array-like, shape = [n_samples, n_classes]
            Returns the log-probabilities of the sample for each class in
            the model, where classes are ordered by arithmetical
            order.

        """
        
        return np.log(self.predict_proba(X))

class DecisionTreeRegressor(BaseDecisionTree, RegressorMixin):
    """Perform regression on dataset with a decision tree.
    
    Example
    -------

    >>> import numpy as np
    >>> from scikits.learn.datasets import load_boston
    >>> from scikits.learn.cross_val import KFold
    >>> from scikits.learn import tree_model
    >>> data = load_boston()
    >>> kf = KFold(len(data.target), 10)
    >>> for train_index, test_index in kf:
    ...     tree = tree_model.DecisionTreeRegressor()
    ...     tree = tree.fit(data.data[train_index], data.target[train_index])
    ...     print np.mean(np.power(tree.predict(data.data[test_index]) - data.target[test_index], 2))
    ... 
    13.3151111111
    7.61977777778
    8.732
    116.885777778
    28.4606666667
    34.4506666667
    13.1571111111
    12.738
    145.118
    29.6389361702

    """
    def __init__(self, criterion='mse', max_depth=10,\
                  min_split=5, F=None, seed=None):       
        BaseDecisionTree.__init__(self, None, 'regression', criterion, \
                                  max_depth, min_split, F, seed)
        

        