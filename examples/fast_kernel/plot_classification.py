"""
============================================================
Comparison of Fast Kernel Machine (Eigenpro) and SVC, part 2
============================================================

Here we train a Fast Kernel Machine (EigenPro) and a Support Vector Classifier (SVC)
on subsets of a synthetic dataset.
Features of this dataset are sampled from two independent Gaussian distributions.
Note that we halt the training for EigenPro in two epochs.
Experimental results
demonstrate that EigenPro achieves high test accuracy, competitive to that of SVC,
while completes training in significant less time.
"""
print(__doc__)


import matplotlib.pyplot as plt
import numpy as np
from time import time

from sklearn.datasets import make_classification
from sklearn.svm import SVC
from sklearn.fast_kernel_classification import FastKernelClassification

rng = np.random.RandomState(1)

centers = np.zeros((2, 50))
centers[0][0] = 2
max_size = 80000
test_size = 10000

# Get data for testing
x, y = make_classification(n_samples=max_size + test_size, n_features=200, n_informative=5, random_state=rng)
x_train = x[:max_size]
y_train = y[:max_size]
x_test = x[max_size:]
y_test = y[max_size:]

fkc_fit_times = []
fkc_pred_times = []
fkc_err = []
svc_fit_times = []
svc_pred_times = []
svc_err = []

train_sizes = [5000, 10000, 20000]

# Fit models to data
for train_size in train_sizes:
    for name, estimator in [("FastKernel", FastKernelClassification(n_epoch=2, bandwidth=10,
                                                                    random_state=1)),
                            ("SupportVector", SVC(C=1, gamma=1. / (2 * 10 * 10), random_state=1))]:
        stime = time()
        estimator.fit(x_train[:train_size], y_train[:train_size])
        fit_t = time() - stime

        stime = time()
        y_pred_test = estimator.predict(x_test)
        pred_t = time() - stime

        err = 100. * np.sum(y_pred_test != y_test) / len(y_test)
        print("Test Error: " + str(err))
        if name == "FastKernel":
            fkc_fit_times.append(fit_t)
            fkc_pred_times.append(pred_t)
            fkc_err.append(err)
        else:
            svc_fit_times.append(fit_t)
            svc_pred_times.append(pred_t)
            svc_err.append(err)
        print("%s Classification with %i training samples in %0.2f seconds." % (name, train_size, fit_t + pred_t))

# Graph time to fit and predict
plt.plot(train_sizes, fkc_fit_times, 'o-', color='r', label='FKC train (fit) time')
plt.plot(train_sizes, svc_fit_times, 'o-', color='g', label='SVC train (fit) time')
plt.plot(train_sizes, fkc_pred_times, 'o--', color='r', label='FKC test (predict) time')
plt.plot(train_sizes, svc_pred_times, 'o--', color='g', label='SVC test (predict) time')

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Train size")
plt.ylabel("Time to fit (seconds)")
plt.title('Classification Fit Time')
plt.legend(loc="best")
plt.xticks(train_sizes)

# Graph accuracy
plt.figure()
plt.plot(train_sizes, fkc_err, 'o-', color='r', label='FastKernel (EigenPro)')
plt.plot(train_sizes, svc_err, 'o-', color='g', label='SupportVector')

plt.xscale("log")
plt.xlabel("Train size")
plt.ylabel("% Error")
plt.title('Accuracy')
plt.legend(loc="best")
plt.show()
