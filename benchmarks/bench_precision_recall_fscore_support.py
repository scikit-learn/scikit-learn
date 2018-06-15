#!/usr/bin/env python
"""
A comparison of multiclass implementation and original implementation for
precision_recall_fscore_support.
"""
from __future__ import print_function

from collections import defaultdict
from time import time

import numpy as np
from numpy import random as nr

from sklearn.metrics.classification import precision_recall_fscore_support, \
    old_precision_recall_fscore_support
import matplotlib.pyplot as plt
from scipy.stats import bernoulli


def compute_bench(samples_range, labels_range):

    it = 0
    results = []

    max_it = len(samples_range) * len(labels_range)
    for n_labels in labels_range:
        res = defaultdict(lambda: [])
        for n_samples in samples_range:
            it += 1
            print('==============================')
            print('Iteration %03d of %03d' % (it, max_it))
            print('==============================')
            print()

            # multilabel case
            y_true = bernoulli.rvs(np.ones((n_samples, n_labels)) / 2,
                                   size=(n_samples, n_labels))
            y_pred = bernoulli.rvs(np.ones((n_samples, n_labels)) / 2,
                                   size=(n_samples, n_labels))

            # multiclass case
            #y_true = nr.randint(0, n_labels, (n_samples,))
            #y_pred = nr.randint(0, n_labels, (n_samples,))

            # binary case
            #y_true = nr.randint(0, 2, (n_samples,))
            #y_pred = nr.randint(0, 2, (n_samples,))
            print("y_true", y_true)

            print('P/R/F/S')
            tstart = time()

            old_precision_recall_fscore_support(y_true, y_pred)

            delta = time() - tstart
            print("Speed: %0.6fs" % delta)
            print()

            res['P/R/F/S'].append(delta)

            print('P/R/F/S with multilabel confusion metrics')

            tstart = time()

            precision_recall_fscore_support(
                y_true, y_pred)

            delta = time() - tstart
            print("Speed: %0.6fs" % delta)
            print()
            print()

            res['P/R/F/S with multilabel confusion metrics'].append(delta)
        results.append(res)

    return results


if __name__ == '__main__':

    samples_range = np.linspace(50, 15000, 20).astype(np.int)
    labels_range = [2, 5, 20, 100, 1000]
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']

    results = compute_bench(samples_range, labels_range)

    plt.figure()
    plt.rcParams["figure.figsize"] = [40, 20]
    for i in range(len(labels_range)):
        plt.plot(samples_range, results[i]['P/R/F/S'],
                 color=colors[i], marker='o',
                 label='Original P/R/F/S, n_label:%d' % labels_range[i])
        plt.plot(samples_range,
                 results[i]['P/R/F/S with multilabel confusion metrics'],
                 color=colors[i], marker='s',
                 label='With multilabel confusion metrics, n_label:%d'
                       % labels_range[i])

    plt.xlim([0, 15000])
    plt.xlabel('Number of samples')
    plt.ylabel('Speed (s)')
    plt.title('precision_recall_fscore_support')
    plt.legend(loc="upper left")
    plt.show()
