# Author: Brian Holt

# License: BSD Style.
"""
================
Random Forest
================

A Random Forest classifier 
    
Implements Random Forests (Breiman 2001)

"""

from __future__ import division
from ..utils import check_random_state
import numpy as np
from scipy import stats
import copy
from ..tree import DecisionTreeClassifier, DecisionTreeRegressor
from ..base import BaseEstimator, ClassifierMixin, RegressorMixin
from ..externals.joblib import Parallel, delayed

__all__ = [
    'RandomForestClassifier',
    'RandomForestRegressor',    
    ]
          
def _train_tree(X, y, i, r, t, random_state):   
    X = np.asanyarray(X, dtype=np.float64, order='C')       
    n_samples, n_features = X.shape    
    y = np.asanyarray(y, dtype=np.float64, order='C') 

    tree = copy.copy(t)

    if r <= 0 or r > 1 :
        raise ValueError("r must be in 0 <= r < 1.\n" +
                         "r is %s" % r)
    n = np.ceil(r * n_samples)
                      
    sample_mask = np.zeros(n_samples, dtype=np.bool)
    permutation = random_state.permutation(n_samples)
    ind_in = permutation[-n:]
    sample_mask[ind_in] = True
    
    X_in = X[sample_mask]
    X_out = X[~sample_mask]
    y_in = y[sample_mask] 
    y_out = y[~sample_mask]
            
    tree.fit(X_in, y_in)
    
    """
    @TODO Compute the out-of-bag error using X_out and y_out
    """
    
    return tree            
            
            
class BaseRandomForest(BaseEstimator):
    '''
    Should not be used directly, use derived classes instead
    '''
    
    def __init__(self, random_state, base_tree,  n_trees, r, n_jobs):
        self.random_state = check_random_state(random_state)       
        self.base_tree = base_tree
        self.n_trees = n_trees
        self.r = r
        self.criterion = base_tree.criterion
        self.max_depth = base_tree.max_depth
        self.min_split = base_tree.min_split
        self.max_features = base_tree.max_features

        if n_jobs <= 0:
            raise ValueError("n_jobs must be >= 0")           
        self.n_jobs = n_jobs
        self.forest = None

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

        self.forest = Parallel(self.n_jobs) \
            (delayed(_train_tree)(X, y, i, self.r, self.base_tree, self.random_state) \
            for i in range(self.n_trees))   
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
            
        if self.forest is None:
            raise Exception('Random Forest not initialized. Perform a fit first')       
        
        C = np.zeros(n_samples)
        for idx, sample in enumerate(X): 
            if self.base_tree.type == 'classification':
                pred_labels = [int(t.predict(sample)) for t in self.forest]
                mode, _ = stats.mode(pred_labels)
                C[idx] = mode.flatten().astype(np.int)
            else:
                pred_labels = [float(t.predict(sample)) for t in self.forest]
                C[idx] = np.mean(pred_labels)           

        return C 


class RandomForestClassifier(BaseRandomForest, ClassifierMixin):
    """Classify a multi-labeled dataset with a random forest.

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

    random_state : integer or array_like, optional
        random_state the random number generator

    n_trees : integer, optional
        the number of trees in the forest
        
    r : float, optional
        the ratio of training samples used per tree 0 < r <= r
    
    n_jobs : integer, optional
        the number of processes to use for parallel computation
    
    #Example
    #-------
    #>>> import numpy as np
    #>>> from scikits.learn.datasets import load_iris
    #>>> from scikits.learn.cross_val import StratifiedKFold
    #>>> from scikits.learn import ensemble
    #>>> data = load_iris()
    #>>> skf = StratifiedKFold(data.target, 10)
    #>>> for train_index, test_index in skf:
    #...     rf = ensemble.RandomForestClassifier(K=3)
    #...     rf.fit(data.data[train_index], data.target[train_index])
    #...     #print np.mean(tree.predict(data.data[test_index]) == data.target[test_index])
    #... 

    
    """     
    def __init__(self, criterion='gini', max_depth=10, min_split=1, 
                 max_features=None, random_state=None, n_trees=10, r=0.7, 
                 n_jobs=1):
        base_tree = DecisionTreeClassifier(criterion=criterion,
            max_depth=max_depth, min_split=min_split,
            max_features=max_features, random_state=random_state)
        BaseRandomForest.__init__(self, random_state, base_tree, n_trees,
                                  r, n_jobs)
           
    def predict_proba(self, X):
        """
        This function does classification or regression on an array of
        test vectors X.

        For a classification model, the probability of each class for each
        sample in X is returned.  

        Parameters
        ----------
        X : array-like, shape = [n_samples, n_features]

        Returns
        -------
        C : array, shape = [n_samples, K]
        """
        
        X = np.atleast_2d(X)
        n_samples, n_features = X.shape    
            
        if self.forest is None:
            raise Exception('Random Forest not initialized. Perform a fit first')       
        
        P = np.zeros((n_samples, self.forest[0].K), dtype=float)
        for t in self.forest:
            P += t.predict_proba(X)
        P /= self.n_trees

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
    
class RandomForestRegressor(BaseRandomForest, RegressorMixin):
    """Perform regression on dataset with a random forest.

    Parameters
    ----------

    criterion : string
        function to measure goodness of split

    max_depth : integer
        maximum depth of the tree  
              
    min_split : integer
        minimum size to split on
        
    F : integer, optional
        if given, then, choose F features

    random_state : integer or array_like, optional
        seed the random number generator

    n_trees : integer, optional
        the number of trees in the forest
        
    r : float, optional
        the ratio of training samples used per tree 0 < r <= r
    
    n_jobs : integer, optional
        the number of processes to use for parallel computation
    
    #Example
    #-------
    #>>> import numpy as np
    #>>> from scikits.learn.datasets import load_boston
    #>>> from scikits.learn.cross_val import KFold
    #>>> from scikits.learn import ensemble
    #>>> data = load_boston()
    #>>> np.random.seed([1]) 
    #>>> perm = np.random.permutation(data.target.size)
    #>>> data.data = data.data[perm]
    #>>> data.target = data.target[perm]    
    #>>> kf = KFold(len(data.target), 2)
    #>>> for train_index, test_index in kf:
    #...     rf = ensemble.RandomForestRegressor(n_jobs=2)
    #...     rf.fit(data.data[train_index], data.target[train_index])
    #...     #print np.mean(np.power(tree.predict(data.data[test_index]) - data.target[test_index], 2)) 
    #... 


   
    """ 
     
    def __init__(self, criterion='mse', max_depth=10, min_split=1,
                 max_features=None, random_state=None, n_trees=10, r=0.7,
                 n_jobs=1):       
        base_tree = DecisionTreeRegressor( criterion=criterion,
            max_depth=max_depth, min_split=min_split,
            max_features=max_features, random_state=random_state)
        BaseRandomForest.__init__(self, random_state, base_tree, n_trees,
                                  r, n_jobs)
