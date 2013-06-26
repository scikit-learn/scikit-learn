import numpy as np
from sklearn.neural_network import SAE
from numpy.testing import assert_almost_equal
from sklearn.datasets import load_digits
X = load_digits().data
X=X.T
sae = SAE()
n_samples=X.shape[1]
n_visible=X.shape[0]
initial_theta=sae._initParams(n_visible)

#This is Unofficial Testing (basic testing): this is run to test if the gradients computed numerically are close enough to those computed analytically

def test_SAEGradient(theta, X,n_visible, n_hidden, lr,sparsityParam, beta, n_samples,n_slice):

   cost,grad=sae._cost(initial_theta, X,n_visible, n_hidden, lr,sparsityParam, beta, n_samples)

   numGrad = computeNumericalGradient(lambda x:  sae._cost(x, X,n_visible, n_hidden, lr,sparsityParam, beta, n_samples)[0], initial_theta,n_slice)
   
   assert_almost_equal(grad[:n_slice], numGrad, decimal=3)
   
def computeNumericalGradient(J, theta,n_slice):
    
    numgrad = np.zeros(n_slice);
    E=np.eye(np.size(theta));
    epsilon=1e-5;
    for i in range(n_slice):
       
        dtheta=E[:,i]*epsilon;
        numgrad[i]=(J(theta+dtheta)-J(theta-dtheta))/epsilon/2.0
    return numgrad


print "Computing analytically and numerically computed gradients ..."
n_slice=20
cost,grad=sae._cost(initial_theta, X,n_visible,n_samples)
numGrad = computeNumericalGradient(lambda x:  sae._cost(x, X,n_visible, n_samples)[0], initial_theta,n_slice)
print 'Analytically Computed Gradient:'
print numGrad
print 'Numerically Computed Gradient:'
print grad[:n_slice]
print "The above matrices should be equal in at least 3 decimal places"




   
