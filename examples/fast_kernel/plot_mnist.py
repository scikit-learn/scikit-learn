"""
============================================================
Comparison of Fast Kernel Machine (EigenPro) and SVC, part 1
============================================================
Here we train a Fast Kernel Machine (EigenPro) and a Support Vector Classifier (SVC)
on subsets of MNIST of various sizes. Note that we halt the training of EigenPro
in two epochs.
We then repeat the same experiments on subsets of MNIST with added label noises.
Specifically, we randomly reset the label (0-9) of 20% samples.
Experimental results on the original MNIST
demonstrate great advantages of EigenPro over SVC in training time.
For MNIST with added label noise, EigenPro shows even larger gain over SVC.
"""
print(__doc__)

import numpy as np
from sklearn.fast_kernel_classification import FastKernelClassification
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from sklearn.datasets import fetch_mldata
from time import time

rng = np.random.RandomState(1)
# Generate sample data from mnist
mnist = fetch_mldata('MNIST original')
mnist.data = mnist.data / 255.

p = np.random.permutation(60000)
x_train = mnist.data[p][:60000]
y_train = np.int32(mnist.target[p][:60000])
x_test = mnist.data[60000:]
y_test = np.int32(mnist.target[60000:])

# randomize 20% of labels
p = np.random.choice(len(y_train), np.int32(len(y_train)*.2), False)
y_train[p] = np.random.choice(10, np.int32(len(y_train)*.2))
p = np.random.choice(len(y_test), np.int32(len(y_test) * .2), False)
y_test[p] = np.random.choice(10, np.int32(len(y_test) * .2))

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
    for name, estimator in [("FastKernel", FastKernelClassification(n_epoch=2, bandwidth=5,random_state=rng)),
                            ("SupportVector", SVC(C=5, gamma=.05))]:
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
        print("%s Classification with %i training samples in %0.2f seconds."%(name, train_size, fit_t+pred_t))


# Graph time
plt.plot(train_sizes, fkc_fit_times, 'o-', color='r', label='FKC train (fit) time')
plt.plot(train_sizes, svc_fit_times, 'o-', color='g', label='SVC train (fit) time')
plt.plot(train_sizes, fkc_pred_times, 'o--', color='r', label='FKC test (predict) time')
plt.plot(train_sizes, svc_pred_times, 'o--', color='g', label='SVC test (predict) time')

plt.xscale("log")
plt.yscale("log")
plt.xlabel("Train size")
plt.ylabel("Time (seconds)")
plt.title('Execution Time (Test size = ' + str(len(x_test)) + ')')
plt.legend(loc="best")

# Graph accuracy
plt.figure()
plt.plot(train_sizes, fkc_err, 'o-', color='r', label='FastKernel')
plt.plot(train_sizes, svc_err, 'o-', color='g', label='SupportVector')

plt.xscale("log")
plt.xlabel("Train size")
plt.ylabel("% Labels guessed incorrectly")
plt.title('Accuracy')
plt.legend(loc="best")
plt.show()
