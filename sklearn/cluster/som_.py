"""
 Self-organizing map 
"""
# Authors: Sebastien Campion <sebastien.campion@inria.fr>
# License: BSD
from __future__ import division
import warnings
from ..base import BaseEstimator
import numpy as np
import math 
import Image 

################################################################################
def save_rgb(neurons,ofile,thumb_size=32):
    '''Function to save map using 3 dim, as RGB'''    
    assert neurons.shape[-1] == 3 
    tsize = (thumb_size,thumb_size) 
    size  = tuple([v * thumb_size for v in neurons.shape[0:2] ])
    im  = Image.new('RGB',size)
    for x in range(neurons.shape[0]):
        for y in range(neurons.shape[1]):
            color = tuple([int(c) for c in neurons[x][y]])
            t = Image.new('RGB',tsize,color)
            im.paste(t,(x*thumb_size,y*thumb_size))            
    im.save(ofile)

################################################################################
class Kohonen2DMap():
    def __init__(self,size,dim,neurons=None):
        #self.log = logging.getLogger('kohonen.map')
        self.dim = dim 
        self.size = size 
        self.neurons = neurons
        if neurons == None  :
            self.neurons = np.random.rand(size,size,dim)#/10
        self.iteration = 0 
        
    def bmu(self,data):
        assert data.shape[0] == self.neurons.shape[-1] 
        data = np.resize(data,self.neurons.shape) 
        dists = np.sum((data-self.neurons)** 2,axis=-1)
        min = dists.argmin()
        #w = np.unravel_index(min,dists.shape)
        return divmod(min,self.size)
        
    def learn(self,datas,nbiter,learning_rate=1,callback=None):
        '''Given an sample of datas, we randomly choose one of them for each 
        iteration.
        A good ratio, nb datas = 2 or 3 x nbiter'''
        self.iteration = 0   
        indices = np.random.random_integers(0,len(datas)-1,nbiter)
        for i in indices: 
            l = nbiter/self.size
            lr = learning_rate * math.exp(-self.iteration/l)
            self._learn_vector(datas[i], nbiter, lr)
            self.iteration += 1 
            if callback != None:
                callback(self,self.iteration)

    def _learn_vector(self, data, nbiter, lr):
        w = self.bmu(data)
        radius = self.radius_of_the_neighbordhood(nbiter)
        for n in self.neurons_in_radius(w,radius):
            nx,ny = n
            wt = self.neurons[nx][ny]
            dr = self.dist(w,n,radius)
            self.neurons[nx][ny] = wt + dr*lr*(data-wt)
            #self.log.debug(('nod',n,'l_rate',lr,'d_radius',dr))
        #self.log.debug(('bmu',w,'iter',self.iteration,'radius',radius))
    
    def dist(self,w,n,radius):
        wx,wy = w
        nx,ny = n
        d = (wx-nx)**2 + (wy-ny)**2
        #offcial paper implementation : return math.exp(-d/2*radius**2)
        return math.exp(-d/radius)
    
    def neurons_in_radius(self,w,radius):
        wi,wj = w 
        r = []
        for i in range(self.neurons.shape[0]):
            for j in range(self.neurons.shape[1]):
                if math.sqrt((i-wi)**2 + (j-wj)**2) < radius:
                    r.append((i,j))
        return r
        
    def radius_of_the_neighbordhood(self,nbiter):
        l = nbiter/self.size
        return self.size * math.exp(-self.iteration/l)
    

################################################################################

class SOM(BaseEstimator):
    """ Self-Organizing Map

    Parameters
    ----------

    data : ndarray
        A M by N array of M observations in N dimensions or a length
        M array of M one-dimensional observations.

    w : int
        Width and height of the square mape as well as the number of
        centroids to generate. If init initialization string is
        'matrix', or if a ndarray is given instead, it is
        interpreted as initial cluster to use instead.

    n_iter : int
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

    def __init__(self, w=16, init='random', n_init=64,learning_rate=1):        
        self.w = w
        self.init = init
        self.n_init = n_init
        self.learning_rate = learning_rate
        self.callback = None

    def fit(self, X, **params):
        """ Compute som"""
        X = np.asanyarray(X)
        self._set_params(**params)

        neurons = None 
        dim = X.shape[-1]
        
        if self.init == 'matrix':
            assert len(self.w.shape) == 3 
            neurons = self.w 
            self.w = neurons.shape[0]

        
        km = Kohonen2DMap(self.w,dim,neurons)
        km.learn(X,self.n_init,self.learning_rate,callback=self.callback)
        self.neurons_ = km.neurons
        self.labels_ = [km.bmu(x) for x in X]
        
        return self
 
