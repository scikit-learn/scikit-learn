"""
A comparison of different optimizers for NCA

Data: UC Irvine's Wine dataset.

"""
import numpy as np

from sklearn.base import clone
from sklearn import metrics
from sklearn.metric_learning import NCATransformer
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import Pipeline
from sklearn.cross_validation import KFold

from sklearn.utils.bench import total_seconds
from sklearn.datasets import fetch_mldata

def get_params(learning_rate=1, loss='kl', max_iter=100, solver='adagrad',
        random_state=0, verbose=2, n_init=10, method='semivectorized',
        threshold=None):
    return {
            "learning_rate" : learning_rate,
            "loss" : loss,
            "max_iter" : max_iter,
            "solver" : solver,
            "random_state" : random_state,
            "verbose" : verbose,
            "n_init" : n_init,
            "method" : method,
            "threshold" : threshold,
        }

if __name__ == '__main__':

    separator = "=" * 30 + "\n"
    wine_ds = fetch_mldata('wine')
    print("Loaded UCI Wine dataset")
    print(separator)

    params = [
            get_params(method='vectorized'),
            get_params(method='semivectorized'),
            get_params(method='vectorized', threshold=0),
            get_params(method='semivectorized', threshold=0),
        ]

    for args in params:
        nca = NCATransformer(**args)
        print(nca)
        nca.fit(wine_ds.data, wine_ds.target)
        print(separator)

