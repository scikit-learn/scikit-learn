import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel

noise=1e-5
signal_var = 8.98576054e+05
length_factor = np.array([5.91326520e+02, 1.32584051e+03])
kernel = ConstantKernel(signal_var, (1e-12, 1e12)) * RBF(length_factor, (1e-12, 1e12)) + WhiteKernel(noise_level=noise)
m = GaussianProcessRegressor(kernel=kernel, alpha=0, n_restarts_optimizer=0, random_state=0)
#Setup Training Data
x_train=np.array([[ 0., 0.],
 [ 1.54919334,-0.77459667],
 [-1.54919334, 0.        ],
 [ 0.         ,-1.54919334],
 [ 0.77459667  ,0.77459667],
 [-0.77459667  ,1.54919334]])
y_train=np.array([[-2.14882017e-10],
 [-4.66975823e+00],
 [ 4.01823986e+00],
 [-1.30303674e+00],
 [-1.35760156e+00],
 [ 3.31215668e+00]])
#fit. We must do this to register teh GPR as "fitted"
m.fit(x_train, y_train)
theta = np.exp(m.kernel_.theta)
x_test = np.array([[-1.93649167,-1.93649167],
 [ 1.93649167, -1.93649167],
 [-1.93649167,  1.93649167],
 [ 1.93649167,  1.93649167]])
#Predict the std_dev in two ways
pred1, std1 = m.predict(x_test, return_std=True)
pred2, cov = m.predict(x_test, return_cov=True)
std2=np.sqrt(np.diagonal(cov))
#Compare
print(std1)
print(std2)