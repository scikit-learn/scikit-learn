""" Incremental Principal Component Analysis
"""

# Author: Kevin Hughes <kevinhughes27@gmail.com>
# License: BSD Style.

import numpy as np
from scipy import linalg as la

from ..base import BaseEstimator, TransformerMixin
from ..utils import array2d, as_float_array
from ..utils.extmath import safe_sparse_dot

class IncrementalPCA(BaseEstimator, TransformerMixin):
    """Incremental principal component analysis (IncrementalPCA)
    
    Linear dimensionality reduction using an online incremental PCA algorithm.
    Components are updated sequentially as new observations are introduced. 
    Each new observation (u) is projected on the eigenspace spanned by
    the current components and the residual vector is used as a new 
    component. The new principal components are then rotated by a rotation 
    matrix (R) whose columns are the eigenvectors of the transformed covariance matrix   
    to yield p + 1 principal components. From those, only the first p are selected.
    
    Based on P.Hall, D. Marshall and R. Martin "Incremental Eigenalysis for 
    Classification" which appeared in British Machine Vision Conference, volume 1,
    pages 286-295, September 1998.
    
    implementation based on the python recipe (recipe-577213-1) 
    published by Micha Kalfon
    
    Parameters
    ----------
    n_components : int
        Number of components to keep.
        Must be set

    copy : bool
        If False, data passed to fit are overwritten

    Attributes
    ----------
    `components_` : array, [n_components, n_features]
        Components with maximum variance.

    `explained_variance_ratio_` : array, [n_components]
        Percentage of variance explained by each of the selected components. \
        k is not set then all components are stored and the sum of explained \
        variances is equal to 1.0
        
    Notes
    -----
    Calling fit(X) multiple times will update the components_ etc.
    
    Examples
    --------

    >>> import numpy as np
    >>> from sklearn.decomposition import IncrementalPCA
    >>> X = np.array([[-1, -1], [-2, -1], [-3, -2], [1, 1], [2, 1], [3, 2]])
    >>> ipca = IncrementalPCA(n_components=2)
    >>> ipca.fit(X)
    IncrementalPCA(copy=True, n_components=2)
    >>> print(ipca.explained_variance_ratio_) # doctest: +ELLIPSIS
    [ 0.99244...  0.00755...]

    See also
    --------
    ProbabilisticPCA
    RandomizedPCA
    KernelPCA
    SparsePCA
    CCIPCA
    """
    def __init__(self, n_components=2, copy=True):
        self.n_components = n_components
        self.copy = copy
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
            self.covariance_ = np.zeros([n_features, n_features], np.float)
            self.components_ = np.zeros([self.n_components,n_features], np.float)
            self.explained_variance_ratio_ = np.zeros([self.n_components], np.float)
        else:
            if n_features != self.components_.shape[1]:
                raise ValueError('The dimensionality of the new data and the existing components_ does not match')   
        
        # incrementally fit the model
        for i in range(0,X.shape[0]):
            self.partial_fit(X[i,:])
        
        # normalize explained_variance_ratio_
        self.explained_variance_ratio_ = (self.explained_variance_ratio_ / self.explained_variance_ratio_.sum())
        
        return self
      
    def partial_fit(self, u):
        """ Updates the mean and components to account for a new vector.

        Parameters
        ----------
        _u : array [1, n_features]
            a single new data sample
        """   

        n = float(self.iteration)
        C = self.covariance_
        V = self.components_
        E = self.explained_variance_ratio_

        # Update covariance matrix and mean vector and centralize input around
        # new mean
        oldmean = self.mean_.copy()
        self.mean_ = (n*self.mean_ + u) / (n + 1.0)
        C = (n*C + np.dot(np.asmatrix(u).T,np.asmatrix(u)) + n*np.dot(np.asmatrix(oldmean).T,np.asmatrix(oldmean)) - (n+1.0)*np.dot(np.asmatrix(self.mean_).T,np.asmatrix(self.mean_))) / (n + 1.0)
        u -= self.mean_    

        # Project new input on current subspace and calculate
        # the normalized residual vector
        g = np.dot(u,V.T)       
        r = u - (np.dot(g,V))
        if np.fabs(r).min() > 1e-9:
            r = (r / np.linalg.norm(r))
        else:
            r = np.zeros_like(r)    

        # Extend the transformation matrix with the residual vector and find
        # the rotation matrix by solving the eigenproblem DR=RE
        V = np.vstack([V, r])
        D = np.dot(V,np.dot(C,V.T))
        (E, R) = np.linalg.eigh(D.T)

        # Sort eigenvalues and eigenvectors from largest to smallest to get the
        # rotation matrix R
        sorter = list(reversed(E.argsort(0)))
        E = E[sorter]
        R = R[:,sorter].T

        # Apply the rotation matrix
        V = np.dot(R,V) 

        # Select only self.n_components largest eigenvectors and values 
        # and update state
        self.iteration += 1
        self.covariance_ = C
        self.components_ = V[0:self.n_components,:]
        self.explained_variance_ratio_ = E[0:self.n_components]

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
        return np.asarray(X_transformed)

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

    `explained_variance_ratio_` : array, [n_components]
        Percentage of variance explained by each of the selected components.

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
    IPCA
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
            self.partial_fit(X[i,:])
        
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
      
    def partial_fit(self, u):
        """ Updates the mean and components to account for a new vector.
        
        Parameters
        ----------
        _u : array [1, n_features]
            a single new data sample
        """
        
        n = float(self.iteration)
        V = self.components_
        
        # amnesic learning params
        if n <= int(self.amnesic):
            w1 = float(n+2-1)/float(n+2)    
            w2 = float(1)/float(n+2)    
        else:
            w1 = float(n+2-self.amnesic)/float(n+2)    
            w2 = float(1+self.amnesic)/float(n+2)

        # update mean
        self.mean_ = w1*self.mean_ + w2*u

        # mean center u        
        u = u - self.mean_

        # update components
        for j in range(0,self.n_components):
            
            if j > n:
                # the component has already been init to a zerovec
                pass
            
            elif j == n:
                # set the component to u 
                V[j,:] = u
            else:       
                # update the components
                V[j,:] = w1*V[j,:] + w2*np.dot(u,V[j,:])*u / la.norm(V[j,:])
                
                normedV = V[j,:] / la.norm(V[j,:])
            
                u = u - np.dot(np.dot(u.T,normedV),normedV)

        self.iteration += 1
        self.components_ = V
            
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
