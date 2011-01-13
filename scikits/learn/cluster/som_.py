"""
 Self-organizing map 
 
 Reference : (to check)
 Kohonen, T.; , "The self-organizing map,"
 Proceedings of the IEEE , vol.78, no.9, pp.1464-1480, Sep 1990 
"""
# Authors: Sebastien Campion <sebastien.campion@inria.fr>
# License: BSD
from __future__ import division
from ..base import BaseEstimator
import numpy as np
import math 

class SelfOrganizingMap(BaseEstimator):
    """ Self-Organizing Map

    Parameters
    ----------

    X : ndarray
        A M by N array of M observations in N dimensions or a length
        M array of M one-dimensional observations.

    size : int
        Width and height of the square map as well as the number of
        centroids to generate. If init initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.

    n_iterations : int
        Number of iterations of the som algrithm to run

    learning_rate : float
        Learning rate

    init : {'random', 'matrix'}
        Method for initialization, defaults to 'random':

        'random': randomly points choosed

        'matrix': interpret the w parameter as a w by M array
         of initial centroids.

    Methods
    -------

    fit(X):
        Compute SOM

    Attributes
    ----------

    neurons_: array, [(x,y), n_features]
        Coordinates of neurons and value

    labels_:
        Labels of each point

    Notes
    ------

    """
    def __init__(self,size=16,init='random',n_iterations=64,learning_rate=1,
                 callback=None ):   
        self.size = size
        self.init = init
        self.n_iterations = n_iterations
        self.learning_rate = learning_rate
        self.callback = callback

    def fit(self, X, **params):
        """Given an sample of X, we randomly choose one of them for each 
        iteration.
        A good ratio, nb X = 2 or 3 x nbiter"""
        X = np.asanyarray(X)
        self._set_params(**params)
        self.dim = X.shape[-1]
        self.neurons_ = None 

        #init neurons_
        if self.init == 'random': 
            self.neurons_ = np.random.rand(self.size,self.size,self.dim)
        elif self.init == 'matrix': 
            assert len(self.size.shape) == 3 
            self.neurons_ = self.size 
            self.size = self.neurons_.shape[0]
    
        #iteration loop 
        self.iteration = 0   
        indices = np.random.random_integers(0,len(X)-1,self.n_iterations)
        for i in indices: 
            l = self.n_iterations/self.size
            lr = self.learning_rate * math.exp(-self.iteration/l)
            self._learn_vector(X[i],lr)
            self.iteration += 1 
            if self.callback != None:
                self.callback(self,self.iteration)
                
        #assign labels
        self.labels_ = [self.bmu(x) for x in X]
        return self

    def _learn_vector(self, vector, lr):
        winner = self.bmu(vector)
        radius = self.radius_of_the_neighbordhood()
        for n in self.neurons_in_radius(winner,radius):
            nx,ny = n
            wt = self.neurons_[nx][ny]
            dr = self.dist(winner,n,radius)
            self.neurons_[nx][ny] = wt + dr*lr*(vector - wt)
 
    def bmu(self,vector):
        """
        best matching unit
        """
        assert vector.shape[0] == self.neurons_.shape[-1] 
        vector = np.resize(vector,self.neurons_.shape) 
        dists = np.sum((vector-self.neurons_)** 2,axis=-1)
        min = dists.argmin()
        #w = np.unravel_index(min,dists.shape)
        return divmod(min,self.size)
    
    def dist(self,w,n,radius):
        wx,wy = w
        nx,ny = n
        d = (wx-nx)**2 + (wy-ny)**2
        #offcial paper implementation : return math.exp(-d/2*radius**2)
        return math.exp(-d/radius)
    
    def neurons_in_radius(self,winner,radius):
        wi,wj = winner 
        r = []
        for i in range(self.neurons_.shape[0]):
            for j in range(self.neurons_.shape[1]):
                if math.sqrt((i-wi)**2 + (j-wj)**2) < radius:
                    r.append((i,j))
        return r
        
    def radius_of_the_neighbordhood(self):
        l = self.n_iterations/self.size
        return self.size * math.exp(-self.iteration/l)
