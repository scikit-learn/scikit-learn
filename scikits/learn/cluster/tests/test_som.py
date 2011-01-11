"""
Testing for SOM 

"""
import Image
import numpy as np
from numpy.testing import assert_equal

from ..som_ import SOM,save_rgb
from .common import generate_clustered_data

n_clusters = 4
X = generate_clustered_data(n_clusters=n_clusters, std=.1)

def test_som():
    np.random.seed(1)
    som = SOM().fit(X,w=4,n_init=32,learning_rate=0.4) 
    labels = som.labels_

    assert_equal(np.unique(labels).shape[0],4)
    assert_equal(np.unique(labels[:20]).shape[0], 1)
    assert_equal(np.unique(labels[20:40]).shape[0], 1)
    assert_equal(np.unique(labels[40:60]).shape[0], 1)
    assert_equal(np.unique(labels[60:]).shape[0], 1)

def test_color_map():    
    train = np.array([[0,0,0],       #black
                      [255,255,255], #white
                      [255,0,0],     #red
                      [0,255,0],     #green
                      [0,0,255],     #blue
                      [255,255,0],   #yellow
                      [0,255,255],   #cyan
                      [255,0,255]    #magenta
                      ])
    w = np.random.rand(16,16,3)*255
    som = SOM(w,n_init=1024,init='matrix',learning_rate=1)
    save_rgb(w,'init.jpg')
    som = SOM(w,n_init=1024,init='matrix',learning_rate=1)
    som.fit(train)
    save_rgb(som.neurons_,'color_map.jpg')


