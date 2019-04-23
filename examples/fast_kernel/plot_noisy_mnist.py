"""
================================================================
Comparison of Fast Kernel Classifier (EigenPro) and SVC on MNIST
with added label noise
================================================================
Here we train a Fast Kernel Machine (EigenPro) and a Support
Vector Classifier (SVC) on subsets of MNIST with added label
noises. Specifically, we randomly reset the label (0-9) of 20%
samples. We halt the training of EigenPro in two epochs.
Experimental results on this dataset demonstrate more than 19 times
speedup of EigenPro over SVC in training time.
EigenPro also shows consistently lower classification error on
test set.
"""

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from time import time

from sklearn.datasets import fetch_mldata
from sklearn.fast_kernel import FKC_EigenPro
from sklearn.svm import SVC

rng = np.random.RandomState(1)
# Generate sample data from mnist
mnist = fetch_mldata('MNIST original')
mnist.data = mnist.data / 255.

p = rng.permutation(60000)
x_train = mnist.data[p][:60000]
y_train = np.int32(mnist.target[p][:60000])
x_test = mnist.data[60000:]
y_test = np.int32(mnist.target[60000:])

# randomize 20% of labels
p = rng.choice(len(y_train), np.int32(len(y_train) * .2), False)
y_train[p] = rng.choice(10, np.int32(len(y_train) * .2))
p = rng.choice(len(y_test), np.int32(len(y_test) * .2), False)
y_test[p] = rng.choice(10, np.int32(len(y_test) * .2))

# Run tests comparing fkc to svc
fkc_fit_times = []
fkc_pred_times = []
fkc_err = []
svc_fit_times = []
svc_pred_times = []
svc_err = []

train_sizes = [500, 1000, 2000]

# Fit models to data
for train_size in train_sizes:
    for name, estimator in [
        ("FastKernel",
         FKC_EigenPro(n_epoch=2, bandwidth=5, random_state=rng)),
            ("SupportVector", SVC(C=5, gamma=1. / (2 * 5 * 5)))]:
        stime = time()
        estimator.fit(x_train[:train_size], y_train[:train_size])
        fit_t = time() - stime

        y_pred_train = estimator.predict(x_train[:train_size])

        stime = time()
        y_pred_test = estimator.predict(x_test)
        pred_t = time() - stime
        train_err = 100. * np.sum(y_pred_train != y_train) / train_size
        err = 100. * np.sum(y_pred_test != y_test) / len(y_test)
        if name == "FastKernel":
            fkc_fit_times.append(fit_t)
            fkc_pred_times.append(pred_t)
            fkc_err.append(err)
        else:
            svc_fit_times.append(fit_t)
            svc_pred_times.append(pred_t)
            svc_err.append(err)
        print("%s Classification with %i training samples in %0.2f seconds. "
              "Train error %.4f" %
              (name, train_size, fit_t + pred_t, train_err))

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
