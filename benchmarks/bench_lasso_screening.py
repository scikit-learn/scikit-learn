"""
Benchmarks of Lasso screening vs other definition of Lasso.

We compare the execution time of Lasso using screening with classical
Lasso and Lasso using precomputed matrices (Gram).

The performances are measured over four differents datasets:
a toy dataset (with a higher number of feature than samples ),
News, MNIST and Leukemia.

We plot the relative time of execution for high and low tolerances.
"""
import gc
import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets

from collections import OrderedDict
from scipy import linalg
from sklearn.datasets.samples_generator import make_regression
from sklearn.linear_model import lasso_path
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.utils.testing import ignore_warnings
from sklearn.exceptions import ConvergenceWarning
from time import time


def load_news():
    data = datasets.fetch_20newsgroups(categories=['comp.graphics',
                                                   'talk.religion.misc'])
    vect = TfidfVectorizer(max_df=0.95, min_df=2, stop_words='english')
    X = vect.fit_transform(data.data)
    y = data.target.astype(np.float)
    y[y == 0] = -1.
    return X, y


def load_leukemia():
    data = datasets.fetch_mldata('leukemia')
    X = data.data
    y = data.target
    X = X.astype(float)
    y = y.astype(float)
    y /= linalg.norm(y)
    return X, y


@ignore_warnings(category=ConvergenceWarning)
def compute_bench(datas, estimators, tol):
    print('Compute benchmark for tol=%.1e' % tol)

    times = np.empty((len(estimators), len(datas)))
    for id_data, (dataset_name, data) in enumerate(datas.items()):
        print('================' + '=' * len(dataset_name))
        print('Dataset %d of %d: %s' % (id_data + 1, len(datas),
                                        dataset_name))
        print('================' + '=' * len(dataset_name))
        X, y = data

        for id_clf, (clf_name, props) in enumerate(estimators.items()):
            gc.collect()
            print("- compute %s" % clf_name)
            tstart = time()
            lasso_path(X, y, eps=1e-2, n_alphas=100, tol=tol,
                       **props)
            t = time() - tstart
            times[id_clf, id_data] = t

    print()
    return times


if __name__ == '__main__':
    # Parameters of the computation
    datas = OrderedDict([
        ('Random p >> n', make_regression(
            n_samples=800, n_features=10000,
            n_informative=100, noise=0.1, random_state=0)),
        ('News', load_news()),
        ('Leukemia', load_leukemia()),
    ])

    estimators = OrderedDict([
        ('Lasso', {
            'precompute': False,
            'screening': 0
        }),
        ('Gram', {
            'precompute': True,
            'screening': 0
        }),
        ('Screening', {
            'precompute': False,
        }),
    ])

    times_table = {tol: compute_bench(datas, estimators, tol=tol)
                   for tol in [1e-2, 1e-8]}

    for tol, times in times_table.items():
        plt.figure(figsize=(10, 6))
        index = np.arange(len(datas))
        bar_width = 0.26
        ptimes = times / np.max(times, axis=0)
        for k, clf_name in enumerate(estimators):
            rects = plt.bar(index + k * bar_width, ptimes[k], bar_width,
                            alpha=0.8, label=clf_name, zorder=10)
            for l, rect in enumerate(rects):
                h = rect.get_height()
                if h > 0.15:
                    h = rect.get_height() - 0.02
                    va = 'top'
                else:
                    h = rect.get_height()
                    va = 'bottom'
                if(times[k, l] > 1):
                    t = '%.1fs' % times[k, l]
                else:
                    t = '%.1fms' % (times[k, l] * 1000)
                plt.text(rect.get_x() + rect.get_width() / 2., h,
                         t, ha='center', va=va, zorder=20)

        plt.ylabel('Relative performances : time / max(time)')
        plt.title('Execution time of Lasso (selection = cyclic, '
                  'tol=%.2e)' % tol)
        plt.xticks(index + bar_width, list(datas))
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.tight_layout()
        plt.subplots_adjust(right=0.865)

    plt.show()
