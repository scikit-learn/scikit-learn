import numpy as N
import pylab as P
from scipy.sandbox.pyem import GM

#------------------------------
# Hyper parameters:
#   - K:    number of clusters
#   - d:    dimension
k   = 3
d   = 2

#-------------------------------------------------------
# Values for weights, mean and (diagonal) variances
#   - the weights are an array of rank 1
#   - mean is expected to be rank 2 with one row for one component
#   - variances are also expteced to be rank 2. For diagonal, one row
#   is one diagonal, for full, the first d rows are the first variance,
#   etc... In this case, the variance matrix should be k*d rows and d 
#   colums
w   = N.array([0.2, 0.45, 0.35])
mu  = N.array([[4.1, 3], [1, 5], [-2, -3]])
va  = N.array([[1, 1.5], [3, 4], [2, 3.5]])

#-----------------------------------------
# First method: directly from parameters:
# Both methods are equivalents.
gm      = GM.fromvalues(w, mu, va)

#-------------------------------------
# Second method to build a GM instance:
gm      = GM(d, k, mode = 'diag')
# The set_params checks that w, mu, and va corresponds to k, d and m
gm.set_param(w, mu, va)

# Once set_params is called, both methods are equivalent. The 2d
# method is useful when using a GM object for learning (where
# the learner class will set the params), whereas the first one
# is useful when there is a need to quickly sample a model
# from existing values, without a need to give the hyper parameters

# Create a Gaussian Mixture from the parameters, and sample
# 1000 items from it (one row = one 2 dimension sample)
data    = gm.sample(1000)

# Plot the samples
P.plot(data[:, 0], data[:, 1], '.')
# Plot the ellipsoids of confidence with a level a 75 %
gm.plot(level = 0.75)
