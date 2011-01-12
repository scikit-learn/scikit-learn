"""
===========================================================
A demo of SelfOrganisingMap with colored neurons 
===========================================================

"""
print __doc__

from time import time
import numpy as np
from scikits.learn.cluster import SelfOrganisingMap
import pylab as pl
import Image 

def map2image(neurons,thumb_size=32):
    '''Function to save som using 3 dim, as RGB'''    
    assert neurons.shape[-1] == 3 
    tsize = (thumb_size,thumb_size) 
    size  = tuple([v * thumb_size for v in neurons.shape[0:2] ])
    im  = Image.new('RGB',size)
    for x in range(neurons.shape[0]):
        for y in range(neurons.shape[1]):
            color = tuple([int(c) for c in neurons[x][y]])
            t = Image.new('RGB',tsize,color)
            im.paste(t,(x*thumb_size,y*thumb_size))            
    return im 

train = np.array([[0,0,0],       #black
                  [255,255,255], #white
                  [255,0,0],     #red
                  [0,255,0],     #green
                  [0,0,255],     #blue
                  [255,255,0],   #yellow
                  [0,255,255],   #cyan
                  [255,0,255]    #magenta
                  ])



init = np.random.rand(16,16,3)*255

pl.subplot(1, 2, 1)
pl.imshow(map2image(init))
pl.title('Initial map')

som = SelfOrganisingMap(init,n_iterations=1024,
                        init='matrix',learning_rate=1)
som.fit(train)

pl.subplot(1, 2, 2)
pl.imshow(map2image(som.neurons_))
pl.title('Organized Map') 
pl.show()
