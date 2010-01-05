from numpy.random import seed

from scipy.sandbox.pyem import GM, GMM, EM
import copy

# To reproduce results, fix the random seed
seed(1)

#+++++++++++++++++++++++++++++
# Meta parameters of the model
#   - k: Number of components
#   - d: dimension of each Gaussian
#   - mode: Mode of covariance matrix: full or diag (string)
#   - nframes: number of frames (frame = one data point = one
#   row of d elements)
k       = 2
d       = 2
mode    = 'diag'
nframes = 1e3

#+++++++++++++++++++++++++++++++++++++++++++
# Create an artificial GM model, samples it
#+++++++++++++++++++++++++++++++++++++++++++
w, mu, va   = GM.gen_param(d, k, mode, spread = 1.5)
gm          = GM.fromvalues(w, mu, va)

# Sample nframes frames  from the model
data    = gm.sample(nframes)

#++++++++++++++++++++++++
# Learn the model with EM
#++++++++++++++++++++++++

# Create a Model from a Gaussian mixture with kmean initialization
lgm = GM(d, k, mode)
gmm = GMM(lgm, 'kmean')

# The actual EM, with likelihood computation. The threshold
# is compared to the (linearly appromixated) derivative of the likelihood
em      = EM()
like    = em.train(data, gmm, maxiter = 30, thresh = 1e-8)

# The computed parameters are in gmm.gm, which is the same than lgm
# (remember, python does not copy most objects by default). You can for example
# plot lgm against gm to compare
