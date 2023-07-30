"""
========================================================================
Benchmark for explicit feature map approximation of polynomial kernels
========================================================================

An example illustrating the approximation of the feature map
of an Homogeneous Polynomial kernel.

.. currentmodule:: sklearn.kernel_approximation

It shows how to use :class:`PolynomialCountSketch` and :class:`Nystroem` to
approximate the feature map of a polynomial kernel for
classification with an SVM on the digits dataset. Results using a linear
SVM in the original space, a linear SVM using the approximate mappings
and a kernelized SVM are compared.

The first plot shows the classification accuracy of Nystroem [2] and
PolynomialCountSketch [1] as the output dimension (n_components) grows.
It also shows the accuracy of a linear SVM and a polynomial kernel SVM
on the same data.

The second plot explores the scalability of PolynomialCountSketch
and Nystroem. For a sufficiently large output dimension,
PolynomialCountSketch should be faster as it is O(n(d+klog k))
while Nystroem is O(n(dk+k^2)). In addition, Nystroem requires
a time-consuming training phase, while training is almost immediate
for PolynomialCountSketch, whose training phase boils down to
initializing some random variables (because is data-independent).

[1] Pham, N., & Pagh, R. (2013, August). Fast and scalable polynomial
kernels via explicit feature maps. In Proceedings of the 19th ACM SIGKDD
international conference on Knowledge discovery and data mining (pp. 239-247)
(http://chbrown.github.io/kdd-2013-usb/kdd/p239.pdf)

[2] Charikar, M., Chen, K., & Farach-Colton, M. (2002, July). Finding frequent
items in data streams. In International Colloquium on Automata, Languages, and
Programming (pp. 693-703). Springer, Berlin, Heidelberg.
(http://www.vldb.org/pvldb/1/1454225.pdf)

"""
# Author: Daniel Lopez-Sanchez <lope@usal.es>
# License: BSD 3 clause

# Load data manipulation functions
# Will use this for timing results
from time import time

# Some common libraries
import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import load_digits
from sklearn.kernel_approximation import Nystroem, PolynomialCountSketch
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline

# Import SVM classifiers and feature map approximation algorithms
from sklearn.svm import SVC, LinearSVC

# Split data in train and test sets
X, y = load_digits()["data"], load_digits()["target"]
X_train, X_test, y_train, y_test = train_test_split(X, y, train_size=0.7)

# Set the range of n_components for our experiments
out_dims = range(20, 400, 20)

# Evaluate Linear SVM
lsvm = LinearSVC().fit(X_train, y_train)
lsvm_score = 100 * lsvm.score(X_test, y_test)

# Evaluate kernelized SVM
ksvm = SVC(kernel="poly", degree=2, gamma=1.0).fit(X_train, y_train)
ksvm_score = 100 * ksvm.score(X_test, y_test)

# Evaluate PolynomialCountSketch + LinearSVM
ps_svm_scores = []
n_runs = 5

# To compensate for the stochasticity of the method, we make n_tets runs
for k in out_dims:
    score_avg = 0
    for _ in range(n_runs):
        ps_svm = Pipeline(
            [
                ("PS", PolynomialCountSketch(degree=2, n_components=k)),
                ("SVM", LinearSVC()),
            ]
        )
        score_avg += ps_svm.fit(X_train, y_train).score(X_test, y_test)
    ps_svm_scores.append(100 * score_avg / n_runs)

# Evaluate Nystroem + LinearSVM
ny_svm_scores = []
n_runs = 5

for k in out_dims:
    score_avg = 0
    for _ in range(n_runs):
        ny_svm = Pipeline(
            [
                (
                    "NY",
                    Nystroem(
                        kernel="poly", gamma=1.0, degree=2, coef0=0, n_components=k
                    ),
                ),
                ("SVM", LinearSVC()),
            ]
        )
        score_avg += ny_svm.fit(X_train, y_train).score(X_test, y_test)
    ny_svm_scores.append(100 * score_avg / n_runs)

# Show results
fig, ax = plt.subplots(figsize=(6, 4))
ax.set_title("Accuracy results")
ax.plot(out_dims, ps_svm_scores, label="PolynomialCountSketch + linear SVM", c="orange")
ax.plot(out_dims, ny_svm_scores, label="Nystroem + linear SVM", c="blue")
ax.plot(
    [out_dims[0], out_dims[-1]],
    [lsvm_score, lsvm_score],
    label="Linear SVM",
    c="black",
    dashes=[2, 2],
)
ax.plot(
    [out_dims[0], out_dims[-1]],
    [ksvm_score, ksvm_score],
    label="Poly-kernel SVM",
    c="red",
    dashes=[2, 2],
)
ax.legend()
ax.set_xlabel("N_components for PolynomialCountSketch and Nystroem")
ax.set_ylabel("Accuracy (%)")
ax.set_xlim([out_dims[0], out_dims[-1]])
fig.tight_layout()

# Now lets evaluate the scalability of PolynomialCountSketch vs Nystroem
# First we generate some fake data with a lot of samples

fakeData = np.random.randn(10000, 100)
fakeDataY = np.random.randint(0, high=10, size=(10000))

out_dims = range(500, 6000, 500)

# Evaluate scalability of PolynomialCountSketch as n_components grows
ps_svm_times = []
for k in out_dims:
    ps = PolynomialCountSketch(degree=2, n_components=k)

    start = time()
    ps.fit_transform(fakeData, None)
    ps_svm_times.append(time() - start)

# Evaluate scalability of Nystroem as n_components grows
# This can take a while due to the inefficient training phase
ny_svm_times = []
for k in out_dims:
    ny = Nystroem(kernel="poly", gamma=1.0, degree=2, coef0=0, n_components=k)

    start = time()
    ny.fit_transform(fakeData, None)
    ny_svm_times.append(time() - start)

# Show results
fig, ax = plt.subplots(figsize=(6, 4))
ax.set_title("Scalability results")
ax.plot(out_dims, ps_svm_times, label="PolynomialCountSketch", c="orange")
ax.plot(out_dims, ny_svm_times, label="Nystroem", c="blue")
ax.legend()
ax.set_xlabel("N_components for PolynomialCountSketch and Nystroem")
ax.set_ylabel("fit_transform time \n(s/10.000 samples)")
ax.set_xlim([out_dims[0], out_dims[-1]])
fig.tight_layout()
plt.show()
