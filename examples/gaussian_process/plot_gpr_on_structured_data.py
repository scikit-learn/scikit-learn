"""
==========================================================================
Gaussian processes on discrete data structures
==========================================================================

This example illustrates the use of Gaussian processes to carry out
regression and classification tasks on data that are not in fixed-length
feature vector form. This is enabled through the use of kernel functions
that can directly operate on discrete structures such as variable-length
sequences, trees, and graphs.
"""
print(__doc__)

import numpy as np
import matplotlib.pyplot as plt
from sklearn.gaussian_process.kernels import Kernel, Hyperparameter
from sklearn.gaussian_process.kernels import StructureOrGenericKernelMixin
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.base import clone


class SequenceKernel(StructureOrGenericKernelMixin, Kernel):
    '''
    a mimimal (but valid) convolutional kernel for sequences of variable length
    '''
    def __init__(self,
                 baseline_similarity=0.5,
                 baseline_similarity_bounds=(1e-5, 1)):
        self.baseline_similarity = baseline_similarity
        self.baseline_similarity_bounds = baseline_similarity_bounds

    @property
    def hyperparameter_baseline_similarity(self):
        return Hyperparameter("baseline_similarity",
                              "numeric",
                              self.baseline_similarity_bounds)

    def _f(self, s1, s2):
        '''
        kernel value between a pair of sequences
        '''
        return sum([1.0 if c1 == c2 else self.baseline_similarity
                   for c1 in s1
                   for c2 in s2])

    def _g(self, s1, s2):
        '''
        kernel derivative between a pair of sequences
        '''
        return sum([0.0 if c1 == c2 else 1.0
                    for c1 in s1
                    for c2 in s2])

    def __call__(self, X, Y=None, eval_gradient=False):
        if Y is None:
            Y = X

        if eval_gradient:
            return (np.array([[self._f(x, y) for y in Y] for x in X]),
                    np.array([[[self._g(x, y)] for y in Y] for x in X]))
        else:
            return np.array([[self._f(x, y) for y in Y] for x in X])

    def diag(self, X):
        return np.array([self._f(x, x) for x in X])

    def is_stationary(self):
        return False

    def clone_with_theta(self, theta):
        cloned = clone(self)
        cloned.theta = theta
        return cloned


kernel = SequenceKernel()

X = ['AGCT', 'AGC', 'AACT', 'TAA', 'AAA', 'GAACA']
y = [1.0, 1.0, 2.0, 2.0, 3.0, 3.0]

'''
Visualize sequence similarity matrix under the kernel
'''

K = kernel(seqs)
D = kernel.diag(seqs)

plt.figure(figsize=(8, 5))
plt.imshow(np.diag(D**-0.5).dot(K).dot(np.diag(D**-0.5)))
plt.gca().set_xticks(np.arange(len(seqs)))
plt.gca().set_xticklabels(seqs)
plt.gca().set_yticks(np.arange(len(seqs)))
plt.gca().set_yticklabels(seqs)
plt.title('Sequence similarity under the kernel')

'''
Regression
'''

training_idx = [0, 1, 3, 4]
gp = GaussianProcessRegressor(kernel)
gp.fit(seqs[training_idx], vals[training_idx])

plt.figure(figsize=(8, 5))
plt.bar(np.arange(len(seqs)), gp.predict(seqs), color='b', label='prediction')
plt.bar(training_idx, vals[training_idx], width=0.2, color='r',
        alpha=0.5, label='training')
plt.gca().set_xticklabels(seqs)
plt.title('Regression on sequences')
plt.legend()

'''
Classification
'''

seqs = np.array(['AGCT', 'CGA', 'TAAC', 'TCG', 'CTTT', 'TGCT'])
# whether there are 'A's in the sequence
clss = np.array([True, True, True, False, False, False])

gp = GaussianProcessClassifier(kernel)
gp.fit(seqs, clss)

seqs_test = ['AAA', 'ATAG', 'CTC', 'CT', 'C']
clss_test = [True, True, False, False, False]

plt.figure(figsize=(8, 5))
plt.scatter(np.arange(len(seqs)), [1.0 if c else -1.0 for c in clss],
            s=100, marker='o', edgecolor='none', facecolor=(1, 0.75, 0),
            label='training')
plt.scatter(len(seqs) + np.arange(len(seqs_test)),
            [1.0 if c else -1.0 for c in clss_test],
            s=100, marker='o', edgecolor='none', facecolor='r', label='truth')
plt.scatter(len(seqs) + np.arange(len(seqs_test)),
            [1.0 if c else -1.0 for c in gp.predict(seqs_test)],
            s=100, marker='x', edgecolor=(0, 1.0, 0.3), linewidth=2,
            label='prediction')
plt.gca().set_xticks(np.arange(len(seqs) + len(seqs_test)))
plt.gca().set_xticklabels(np.concatenate((seqs, seqs_test)))
plt.gca().set_yticks([-1, 1])
plt.gca().set_yticklabels([False, True])
plt.title('Classification on sequences')
plt.legend()

plt.show()
