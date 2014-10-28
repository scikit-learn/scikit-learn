"""
===========================
Multi-class classification
===========================

An example show how to use the different encoding methods in scikit-learn for
multi-class problem, for example, a hand-written digits recognition. And a
comparison is conducted among different encoding methods in scikit-learn:
(1) Using the label_binarizer in SVC; (2) OneVsOneClassifier;
(3) OneVsRestClassifier; (4) OutputCodeClassifer with two strategy, 'random'
and 'iter_hamming'.

"""
print(__doc__)

# Author: Qichao Que <que@cse.ohio-state.edu>
# License: BSD 3 clause

# Import timing, plotting and scientific python libarary.
from matplotlib import pyplot as plt
import numpy as np
import time

# Import datasets, classifiers, and performance metrics
from sklearn.svm import SVC
from sklearn import datasets, metrics
from sklearn import multiclass
from sklearn.metrics.pairwise import pairwise_distances

# Load the digits dataset
digits = datasets.load_digits()

# Split the dataset in to training and testing data
train_X = digits.data[:1000]
train_y = digits.target[:1000]
test_X = digits.data[1000:]
test_y = digits.target[1000:]

# Using label_binarizer built-in SVC
c = SVC(gamma=0.001)
s = time.time()
c.fit(train_X, train_y)
pred = c.predict(test_X)
e = time.time()
time_svc = e - s
error_svc = np.sum((test_y != pred).astype(np.int)) / test_y.shape[0]
print('Error for SVC: %.3f%%' % error_svc * 100)

# Using OneVsRestClassifier
c1 = multiclass.OneVsRestClassifier(c)
s = time.time()
c1.fit(train_X, train_y)
pred = c1.predict(test_X)
e = time.time()
time_ovr = e - s
error_ovr = np.sum((test_y != pred).astype(np.int)) / test_y.shape[0]
print('Error for SVC using OneVsRest Code: %.3f%%' % error_ovr * 100)

# Using OneVsOneClassifier
c1 = multiclass.OneVsOneClassifier(c)
s = time.time()
c1.fit(train_X, train_y)
pred = c1.predict(test_X)
e = time.time()
time_ovo = e - s
error_ovo = np.sum((test_y != pred).astype(np.int)) / test_y.shape[0]
print('Error for SVC using OneVsOne Code: %.3f%%' % error_ovo * 100)

# constants for OutputCodeClassifier
max_iter = 50
repeat = 20
code_size = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

# Using OutputCodeClassifier with strategy 'iter_hamming'
error_ecoc_iter_hamming = np.zeros((len(code_size), repeat))
random_state = np.random.RandomState(0)
time_ecoc_iter_hamming = np.zeros((len(code_size), repeat))
for j, cs in enumerate(code_size):
  for i in range(repeat):
    c1 = multiclass.OutputCodeClassifier(c, max_iter=max_iter, code_size=cs,
                                         strategy="iter_hamming",
                                         random_state=random_state)
    s = time.time()
    c1.fit(train_X, train_y)
    pred = c1.predict(test_X)
    e = time.time()
    time_ecoc_iter_hamming[j, i] = e - s
    error_ecoc_iter_hamming[j, i] = np.sum((test_y != pred).astype(np.int)) / test_y.shape[0]
mean_error_ecoc_iter_hamming = np.mean(error_ecoc_iter_hamming, axis=1)
print('Error for SVC using Output Code with iter_hamming strategy: ')
print(mean_error_ecoc_iter_hamming)

# Using OutputCodeClassifier with strategy 'random'
error_ecoc_random = np.zeros((len(code_size), repeat))
random_state = np.random.RandomState(0)
time_ecoc_random = np.zeros((len(code_size), repeat))
for j, cs in enumerate(code_size):
  for i in range(repeat):
    c1 = multiclass.OutputCodeClassifier(c, max_iter=max_iter, code_size=cs,
                                         strategy="random",
                                         random_state=random_state)
    s = time.time()
    c1.fit(train_X, train_y)
    pred = c1.predict(test_X)
    e = time.time()
    time_ecoc_random[j, i] = e-s
    error_ecoc_random[j, i] = np.sum((test_y!=pred).astype(np.int))/test_y.shape[0]
mean_error_ecoc_random = np.mean(error_ecoc_random, axis=1)
print('Error for SVC using Output Code with random strategy:')
print(mean_error_ecoc_random)

plt.figure(1)
# Plot Result, the classification error using different encoding strategy.
plt.subplot(121)
plt.plot(1, error_svc, 'b', marker='+', label='SVC')
plt.plot(1, error_ovo, 'r', marker='o', label='OneVsOne')
plt.plot(1, error_ovr, 'g', marker='*', label='OneVsRest')
plt.plot(code_size, np.mean(error_ecoc_iter_hamming, axis=1), '-c', marker='x',
         label='OC-iter_hamming')
plt.plot(code_size, np.mean(error_ecoc_random, axis=1), '-k', marker='s',
         label='OC-random')
plt.legend()
plt.xlabel('code_size')
plt.ylabel('Classification Error')

# Plot the timing results as well.
plt.subplot(122)
plt.plot(1, time_svc, 'b', marker='+', label='SVC')
plt.plot(1, time_ovo, 'r', marker='o', label='OneVsOne')
plt.plot(1, time_ovr, 'g', marker='*', label='OneVsRest')
plt.plot(code_size, np.mean(time_ecoc_iter_hamming, axis=1), '-c', marker='x',
         label='iter_hamming')
plt.plot(code_size, np.mean(time_ecoc_random, axis=1), '-k', marker='s',
         label='random')
plt.xlabel('code_size')
plt.ylabel('Running time')
plt.legend()

plt.show()
