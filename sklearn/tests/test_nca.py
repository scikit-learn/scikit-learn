import numpy as np
from sklearn.metric_learning import NCA
from sklearn.utils.testing import ignore_warnings
from sklearn.utils.testing import assert_array_almost_equal

N = 300
aux = (np.concatenate([0.5*np.ones((N/2, 1)),
                       np.zeros((N/2, 1)), 1.1*np.ones((N/2, 1))], axis=1))
X = np.concatenate([np.random.rand(N/2, 3),
                    np.random.rand(N/2, 3) + aux])

y = np.concatenate([np.concatenate([np.ones((N/2, 1)), np.zeros((N/2, 1))]),
                    np.concatenate([np.zeros((N/2, 1)), np.ones((N/2, 1))])],
                   axis=1)
X = X.T
y = y[:, 0]
A = np.array([[1, 0, 0], [0, 1, 0]])

def test_NCA():
	# Training
	nca = NCA.NCA(metric=A, method='BFGS', objective='kl-divergence', options={'maxiter': 10, 'disp': True})
	print nca.score(X, y)
	nca.fit(X, y)
	print nca.score(X, y)