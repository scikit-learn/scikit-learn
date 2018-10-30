"""
=======================================================
Comparison of Fast Kernel Machine (Eigenpro) and SVC on
synthetic dataset
=======================================================
Here we train a Fast Kernel Machine (EigenPro) and a
Support Vector Classifier (SVC) on subsets of a synthetic
dataset. Features of this dataset are sampled from two
independent Gaussian distributions. We halt the training
for EigenPro in 3 epochs. Experimental results
demonstrate that EigenPro achieves high test accuracy,
competitive to that of SVC, while completes training in
significant less time (8 times speedup).
"""
print(__doc__)

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from time import time

from sklearn.datasets import make_classification
from sklearn.fast_kernel import FKC_EigenPro
from sklearn.svm import SVC

rng = np.random.RandomState(1)

centers = np.zeros((2, 50))
centers[0][0] = 2
max_size = 10000
test_size = 10000

# Get data for testing
x, y = make_classification(n_samples=max_size+test_size,
                           n_features=200,
                           n_informative=10,
                           random_state=rng)
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

train_sizes = [2000, 5000, 10000]

for train_size in train_sizes:
    for name, estimator in [
        ("FastKernel", FKC_EigenPro(n_epoch=3, bandwidth=10,
                                    random_state=rng)),
            ("SupportVector", SVC(C=5, gamma=1./(2 * 10 * 10)))]:
        stime = time()
        estimator.fit(x_train[:train_size], y_train[:train_size])
        fit_t = time() - stime

        stime = time()
        y_pred_test = estimator.predict(x_test)
        pred_t = time() - stime

        err = 100. * np.sum(y_pred_test != y_test) / len(y_test)
        if name == "FastKernel":
            fkc_fit_times.append(fit_t)
            fkc_pred_times.append(pred_t)
            fkc_err.append(err)
        else:
            svc_fit_times.append(fit_t)
            svc_pred_times.append(pred_t)
            svc_err.append(err)
        print("%s Classification with %i training samples in %0.2f seconds." %
              (name, train_size, fit_t + pred_t))

# set up grid for figures
fig = plt.figure(num=None, figsize=(6, 4), dpi=160)
ax = plt.subplot2grid((2, 2), (0, 0), rowspan=2)
train_size_labels = [str(s) for s in train_sizes]

# Graph fit(train) time
ax.plot(train_sizes, svc_fit_times, 'o--', color='g', label='SVC')
ax.plot(train_sizes, fkc_fit_times, 'o-', color='r', label='FKC (EigenPro)')
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ax.set_xlabel('train size')
ax.set_ylabel('time (seconds)')

ax.legend()
ax.set_title('Train set')
ax.set_xticks(train_sizes)
ax.set_xticklabels(train_size_labels)
ax.set_xticks([], minor=True)
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

# Graph prediction(test) time
ax = plt.subplot2grid((2, 2), (0, 1), rowspan=1)
ax.plot(train_sizes, fkc_pred_times, 'o-', color='r')
ax.plot(train_sizes, svc_pred_times, 'o--', color='g')
ax.set_xscale('log')
ax.set_yscale('log', nonposy='clip')
ax.set_ylabel('time (seconds)')
ax.set_title('Test set')
ax.set_xticks([])
ax.set_xticks([], minor=True)

# Graph training error
ax = plt.subplot2grid((2, 2), (1, 1), rowspan=1)
ax.plot(train_sizes, fkc_err, 'o-', color='r')
ax.plot(train_sizes, svc_err, 'o-', color='g')
ax.set_xscale('log')
ax.set_xticks(train_sizes)
ax.set_xticklabels(train_size_labels)
ax.set_xticks([], minor=True)
ax.set_xlabel('train size')
ax.set_ylabel('classification error %')
plt.tight_layout()
plt.show()
