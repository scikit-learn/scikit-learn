import time

import numpy as np

from sklearn.svm import SVR
from sklearn.grid_search import GridSearchCV
from sklearn.learning_curve import learning_curve
from sklearn.kernel_ridge import KernelRidge
import matplotlib.pyplot as plt

np.random.seed(0)

###############################################################################
# Generate sample data
X = 5 * np.random.rand(10000, 1)
y = np.sin(X).ravel()

# Add noise to targets
y[::5] += 3 * (0.5 - np.random.rand(X.shape[0]/5))

X_plot = np.linspace(0, 5, 100000)[:, None]

###############################################################################
# Fit regression model
svr = GridSearchCV(SVR(kernel='rbf', gamma=1e-1), cv=5,
				   param_grid={"C": [1e0, 1e1, 1e2, 1e3],
							   "gamma": np.logspace(-2, 2, 5)
							   })

kr = GridSearchCV(KernelRidge(kernel='rbf', gamma=1e-1), cv=5,
				  param_grid={"alpha": [1e0, 1e-1, 1e-2, 1e-3],
				  			  "gamma": np.logspace(-2, 2, 5)
				  			  })

t0 = time.time()
svr.fit(X[:100], y[:100])
print "SVR complexity and bandwidth selected and model fitted in %.3f s" \
	% (time.time() - t0)

t0 = time.time()
kr.fit(X[:100], y[:100])
print "KR complexity and bandwidth selected and model fitted in %.3f s" \
	% (time.time() - t0)

t0 = time.time()
y_svr = svr.predict(X_plot)
print "SVR prediction for %d inputs in %.3f s" \
	% (X_plot.shape[0], time.time() - t0)

t0 = time.time()
y_kr = kr.predict(X_plot)
print "KR prediction for %d inputs in %.3f s" \
	% (X_plot.shape[0], time.time() - t0)


###############################################################################
# look at the results
plt.scatter(X[:100], y[:100], c='k', label='data')
plt.hold('on')
plt.plot(X_plot, y_svr, c='r', label='SVR')
plt.plot(X_plot, y_kr, c='g', label='KR')
plt.xlabel('data')
plt.ylabel('target')
plt.title('Support Vector Regression')
plt.legend()

# Visualize training and prediction time
plt.figure()

# Generate sample data
X = 5 * np.random.rand(10000, 1)
y = np.sin(X).ravel()
y[::5] += 3 * (0.5 - np.random.rand(X.shape[0]/5))
sizes = np.logspace(1, 4, 7)
for name, estimator in {"kr": KernelRidge(kernel='rbf', alpha=1e-1, gamma=1e+1),
             		    "svr": SVR(kernel='rbf', C=1e1, gamma=1e+1)}.items():
	train_time = []
	test_time = []
	for train_test_size in sizes:
		t0 = time.time()
		estimator.fit(X[:train_test_size], y[:train_test_size])
		train_time.append(time.time() - t0)

		t0 = time.time()
		estimator.predict(X_plot[:train_test_size])
		test_time.append(time.time() - t0)

	plt.plot(sizes, train_time, 'o-', color="r" if name == "svr" else "g",
 		 label="%s (train)" % name)
	plt.plot(sizes, test_time, 'o--', color="r" if name == "svr" else "g",
		 label="%s (test)" % name)

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Train/Test size")
plt.ylabel("Time (seconds)")
plt.title('Execution Time')
plt.legend(loc="best")

# Visualize learning curves
plt.figure()

svr = SVR(kernel='rbf', C=1e1, gamma=1e-1)
kr = KernelRidge(kernel='rbf', alpha=1e-1, gamma=1e-1)
train_sizes, train_scores_svr, test_scores_svr = \
	learning_curve(svr, X[:100], y[:100], train_sizes=np.linspace(0.1, 1, 10),
 				   scoring="mean_squared_error", cv=10)
train_sizes_abs, train_scores_kr, test_scores_kr = \
 	learning_curve(kr, X[:100], y[:100], train_sizes=np.linspace(0.1, 1, 10),
 				   scoring="mean_squared_error", cv=10)

plt.plot(train_sizes, test_scores_svr.mean(1), 'o-', color="r",
 		 label="SVR")
plt.plot(train_sizes, test_scores_kr.mean(1), 'o-', color="g",
		 label="Kernel Regression")
plt.xlabel("Train size")
plt.ylabel("Mean Squared Error")
plt.title('Learning curves')
plt.legend(loc="best")

plt.show()
