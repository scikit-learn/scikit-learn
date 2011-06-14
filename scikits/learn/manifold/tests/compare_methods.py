import numpy
import pylab
from time import time

from scikits.learn.manifold import locally_linear_embedding
from scikits.learn.datasets.samples_generator import s_curve, swiss_roll

X,c = s_curve(1000)
#X,c = swiss_roll(1000)

k = 8
d_out = 2

methods = ['standard', 'hessian', 'modified']

pylab.figure(figsize=(10,8))
for i,method in enumerate(methods):
    t0 = time()
    Y,err  = locally_linear_embedding(X,k,d_out,
                                      eigen_solver='arpack',
                                      method=method)
    t1 = time()
    print "%s: %.2g sec" % (methods[i],t1-t0)
    print ' err = %.2e' % err
    
    pylab.subplot(222+i)
    pylab.scatter(Y[:,0],Y[:,1],c=c)
    pylab.title(methods[i])

pylab.show()
