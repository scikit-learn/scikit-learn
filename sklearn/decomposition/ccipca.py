""" Candid Covariance-free Incremental Principal Component Analysis
"""

# Author: Kevin Hughes <kevinhughes27@gmail.com>
# License: BSD Style.

import numpy as np
from scipy import linalg as la

from ..base import BaseEstimator, TransformerMixin
from ..utils import array2d, as_float_array
from ..utils.extmath import safe_sparse_dot

class CCIPCA(BaseEstimator, TransformerMixin):
    """Candid covariance-free incremental principal component analysis (CCIPCA)

    Linear dimensionality reduction using an online incremental PCA algorithm.
    CCIPCA computes the principal components incrementally without
    estimating the covariance matrix. This algorithm was designed for high
    dimensional data and converges quickly. 

    This implementation only works for dense arrays. However it should scale
    well to large data.
    
    Time Complexity: per iteration 3 dot products and 2 additions over 'n' where 'n' 
                     is the number of features (n_features).

    Implementation of:
        author={Juyang Weng and Yilu Zhang and Wey-Shiuan Hwang}
        journal={Pattern Analysis and Machine Intelligence, IEEE Transactions on}
        title={Candid covariance-free incremental principal component analysis}
        year={2003}
        month={aug}
        volume={25}
        number={8}
        pages={1034-1040}

    Parameters
    ----------
    n_components : int
        Number of components to keep.
        Must be set
    
    amnesia : float
        A parameter that weights the present more strongly than the
        past. amnesia=1 makes the present count the same as anything else.

    copy : bool
        If False, data passed to fit are overwritten

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Components.

    Notes
    -----
    Calling fit(X) multiple times will update the components_ etc.

    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.decomposition import CCIPCA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> ccipca = CCIPCA(n_components=2)
    >>> ccipca.fit(X)
    CCIPCA(amnesic=2.0, copy=True, n_components=2)
    >>> print(ccipca.explained_variance_ratio_)
    [ 0.97074203  0.02925797]
    
    See also
    --------
    PCA
    ProbabilisticPCA
    RandomizedPCA
    KernelPCA
    SparsePCA
    """
    
    def __init__(self, n_components=2, amnesic=2.0, copy=True):
        self.n_components = n_components
        if self.n_components < 2:
            raise ValueError ("must specifiy n_components for CCIPCA")
            
        self.copy = copy
        self.amnesic = amnesic
        self.iteration = 0

    def fit(self, X, y=None, **params):
        """Fit the model with X.

        Parameters
        ----------
        X: array-like, shape (n_samples, n_features)
            Training data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        self : object
            Returns the instance itself.
            
        Notes
        -----
        Calling multiple times will update the components
        """
        
        X = array2d(X)
        n_samples, n_features = X.shape 
        X = as_float_array(X, copy=self.copy)
        
        # init
        if self.iteration == 0:  
            self.mean_ = np.zeros([n_features], np.float)
            self.components_ = np.zeros([self.n_components,n_features], np.float)
        else:
            if n_features != self.components_.shape[1]:
                raise ValueError('The dimensionality of the new data and the existing components_ does not match')   
        
        # incrementally fit the model
        for i in range(0,X.shape[0]):
            self.inc_fit(X[i,:])
        
        # update explained_variance_ratio_
        self.explained_variance_ratio_ = np.sqrt(np.sum(self.components_**2,axis=1))
        
        # sort by explained_variance_ratio_
        idx = np.argsort(-self.explained_variance_ratio_)
        self.explained_variance_ratio_ = self.explained_variance_ratio_[idx]
        self.components_ = self.components_[idx,:]
        
        # re-normalize
        self.explained_variance_ratio_ = (self.explained_variance_ratio_ / self.explained_variance_ratio_.sum())
            
        for r in range(0,self.components_.shape[0]):
            self.components_[r,:] /= np.sqrt(np.dot(self.components_[r,:],self.components_[r,:]))
        
        return self
      
    def inc_fit(self, u):
        """ Updates the mean and components to account for a new vector.
        
        Parameters
        ----------
        _u : array [1, n_features]
            a single new data sample
        """
        
        # amnesic learning params
        if self.iteration <= int(self.amnesic):
            w1 = float(self.iteration+2-1)/float(self.iteration+2)    
            w2 = float(1)/float(self.iteration+2)    
        else:
            w1 = float(self.iteration+2-self.amnesic)/float(self.iteration+2)    
            w2 = float(1+self.amnesic)/float(self.iteration+2)

        # update mean
        self.mean_ = w1*self.mean_ + w2*u

        # mean center u        
        u = u - self.mean_

        # update components
        for j in range(0,self.n_components):
            
            if j > self.iteration:
                # the component has already been init to a zerovec
                pass
            
            elif j == self.iteration:
                # set the component to u 
                self.components_[j,:] = u
            else:       
                # update the components
                self.components_[j,:] = w1*self.components_[j,:] + w2*np.dot(u,self.components_[j,:])*u / la.norm(self.components_[j,:])
                
                normedV = self.components_[j,:] / la.norm(self.components_[j,:])
            
                u = u - np.dot(np.dot(u.T,normedV),normedV)

        self.iteration += 1
            
        return
    
    def transform(self, X):
        """Apply the dimensionality reduction on X.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            New data, where n_samples in the number of samples
            and n_features is the number of features.

        Returns
        -------
        X_new : array-like, shape (n_samples, n_components)

        """
        X = array2d(X)
        X_transformed = X - self.mean_
        X_transformed = np.dot(X_transformed, self.components_.T)
        return X_transformed

    def inverse_transform(self, X):
        """Transform data back to its original space, i.e.,
        return an input X_original whose transform would be X

        Parameters
        ----------
        X : array-like, shape (n_samples, n_components)
            New data, where n_samples in the number of samples
            and n_components is the number of components.

        Returns
        -------
        X_original array-like, shape (n_samples, n_features)
        """
        return np.dot(X, self.components_) + self.mean_
