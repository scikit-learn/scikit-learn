import pylab as pl
import numpy as np
import scipy
from scikits.learn.datasets.samples_generator import multivariate_normal_from_latent_variables


# 
latent_coefs = np.array([[1]*5 + [0]*15,
                                 [0]*5  + [1]*5 + [0]*10,
                                 [0]*10 + [1]*5 + [0]*5,
                                 [0]*15 + [1]*5,
                                 [0]*5  + [1]*10 + [0]*5])
print latent_coefs
M, latent = multivariate_normal_from_latent_variables(latent_coefs)
pl.matshow(np.corrcoef(M, rowvar=0)); pl.colorbar()

## Param
X = M[:,0:10]
Y = M[:,11:19]
n_components = 2

run /home/duchesnay/git/scikit-learn/scikits/learn/pls.py

pls = PLS(n_components=n_components)
pls.fit(X,Y)

C = np.dot(X.T,Y)
U, s, Vh = linalg.svd(C)

[scipy.correlate(U[:,k],pls.u_[:,k]) for k in 1:n_components]

cor_u = [scipy.correlate(U[:,k],pls.u_[:,k])[0] for k in xrange(n_components)]
cor_v = [scipy.correlate(Vh.T[:,k],pls.v_[:,k])[0] for k in xrange(n_components)]
np.all(np.abs(cor_u)>.99)
np.all(np.abs(cor_v)>.99)

