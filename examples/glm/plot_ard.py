"""
==================================================
Automatic Relevance Determination Regression
==================================================
"""

from scikits.learn.glm import ARDRegression
import numpy as np
import pylab as pl
import numpy.random as nr
import scipy.stats as st


################################################################################
# Generating simulated data with Gaussian weigthts

### Parameters of the example
nr.seed(0)
n_samples = 50
n_features = 100
### Create gaussian data
X = nr.randn(n_samples, n_features)
### Create weigts with a precision lambda_ of 4.
lambda_ = 4.
w = np.zeros(n_features)
### Only keep 10 weights of interest
relevant_features = nr.randint(0,n_features,10)
for i in relevant_features:
  w[i] = st.norm.rvs(loc = 0, scale = 1./np.sqrt(lambda_))
### Create noite with a precision alpha of 50.
alpha_ = 50.
noise =  st.norm.rvs(loc = 0, scale = 1./np.sqrt(alpha_), size = n_samples)
### Create the target
Y = np.dot(X, w) + noise


################################################################################
### Fit the ARD Regression
clf = ARDRegression(compute_score = True)
clf.fit(X, Y)



################################################################################
### Plot the true weights, the estimated weights and the histogram of the
### weights

pl.figure()
axe = pl.axes([0.1,0.6,0.8,0.325])
axe.set_title("ARD - Weights of the model")
axe.plot(clf.coef_, 'b-', label="Estimate")
axe.plot(w, 'g-', label="Ground truth")
axe.set_xlabel("Features")
axe.set_ylabel("Values of the weights")
axe.legend(loc=1)

axe = pl.axes([0.1,0.1,0.45,0.325])
axe.set_title("Histogram of the weights")
axe.hist(clf.coef_, bins=n_features, log=True)
axe.plot(clf.coef_[relevant_features],5*np.ones(len(relevant_features)),'ro',
label="Relevant features")
axe.set_ylabel("Features")
axe.set_xlabel("Values of the weights")
axe.legend(loc=1)

axe = pl.axes([0.65,0.1,0.3,0.325])
axe.set_title("Objective function")
axe.plot(clf.all_score_)
axe.set_ylabel("Score")
axe.set_xlabel("Iterations")
pl.show()

