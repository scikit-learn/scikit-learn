from urllib.request import urlretrieve
import os
from gzip import GzipFile
from time import time
import argparse

import numpy as np
import pandas as pd
from joblib import Memory
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.ensemble import FastGradientBoostingClassifier
from sklearn.ensemble import FastGradientBoostingRegressor
from sklearn.datasets import make_classification
from sklearn.datasets import make_regression
from sklearn._fast_gradient_boosting.utils import get_lightgbm_estimator


parser = argparse.ArgumentParser()
parser.add_argument('--n-leaf-nodes', type=int, default=31)
parser.add_argument('--n-trees', type=int, default=10)
parser.add_argument('--lightgbm', action="store_true", default=False,
                    help='also plot lightgbm')
parser.add_argument('--learning-rate', type=float, default=.1)
parser.add_argument('--problem', type=str, default='classification',
                    choices=['classification', 'regression'])
parser.add_argument('--n-classes', type=int, default=2)
parser.add_argument('--n-samples-max', type=int, default=int(1e6))
parser.add_argument('--n-features', type=int, default=20)
parser.add_argument('--max-bins', type=int, default=255)
args = parser.parse_args()

n_leaf_nodes = args.n_leaf_nodes
n_trees = args.n_trees
lr = args.learning_rate
max_bins = args.max_bins

def get_estimator_and_data():
    if args.problem == 'classification':
        X, y = make_classification(args.n_samples_max,
                                   n_features=args.n_features,
                                   n_classes=args.n_classes,
                                   n_clusters_per_class=1,
                                   random_state=0)
        return X, y, FastGradientBoostingClassifier
    elif args.problem == 'regression':
        X, y = make_regression(args.n_samples_max,
                               n_features=args.n_features, random_state=0)
        return X, y, FastGradientBoostingRegressor


X, y, Estimator = get_estimator_and_data()
X_train_, X_test_, y_train_, y_test_ = train_test_split(X, y, random_state=0)


def one_run(n_samples):
    X_train = X_train_[:n_samples]
    X_test = X_test_[:n_samples]
    y_train = y_train_[:n_samples]
    y_test = y_test_[:n_samples]

    print("Fitting a sklearn model...")
    tic = time()
    est = Estimator(learning_rate=lr,
                    n_estimators=n_trees,
                    max_bins=max_bins,
                    max_leaf_nodes=n_leaf_nodes,
                    n_iter_no_change=None,
                    random_state=0,
                    verbose=0)
    est.fit(X_train, y_train)
    sklearn_fit_duration = time() - tic
    tic = time()
    sklearn_score = est.score(X_test, y_test)
    sklearn_score_duration = time() - tic
    print("score: {:.4f}".format(sklearn_score))
    print("fit duration: {:.3f}s,".format(sklearn_fit_duration))
    print("score duration: {:.3f}s,".format(sklearn_score_duration))

    if args.lightgbm:
        print("Fitting a LightGBM model...")
        # get_lightgbm does not accept loss='auto'
        if args.problem == 'classification':
            loss = 'binary_crossentropy' if args.n_classes == 2 else \
                'categorical_crossentropy'
            est.set_params(loss=loss)
        lightgbm_est = get_lightgbm_estimator(est)

        tic = time()
        lightgbm_est.fit(X_train, y_train)
        lightgbm_fit_duration = time() - tic
        tic = time()
        lightgbm_score = lightgbm_est.score(X_test, y_test)
        lightgbm_score_duration = time() - tic
        print("score: {:.4f}".format(lightgbm_score))
        print("fit duration: {:.3f}s,".format(lightgbm_fit_duration))
        print("score duration: {:.3f}s,".format(lightgbm_score_duration))

        return (sklearn_score, sklearn_fit_duration, sklearn_score_duration,
                lightgbm_score, lightgbm_fit_duration,
                lightgbm_score_duration)

    return (sklearn_score, sklearn_fit_duration, sklearn_score_duration,
            None, None, None)

n_samples_list = [1000, 10000, 100000, 500000, 1000000, 5000000, 10000000]
n_samples_list = [n_samples for n_samples in n_samples_list
                  if n_samples <= args.n_samples_max]

sklearn_scores = []
sklearn_fit_durations = []
sklearn_score_durations = []
lightgbm_scores = []
lightgbm_fit_durations = []
lightgbm_score_durations = []

for n_samples in n_samples_list:
    (sklearn_score,
     sklearn_fit_duration,
     sklearn_score_duration,
     lightgbm_score,
     lightgbm_fit_duration,
     lightgbm_score_duration) = one_run(n_samples)

    sklearn_scores.append(sklearn_score)
    sklearn_fit_durations.append(sklearn_fit_duration)
    sklearn_score_durations.append(sklearn_score_duration)
    lightgbm_scores.append(lightgbm_score)
    lightgbm_fit_durations.append(lightgbm_fit_duration)
    lightgbm_score_durations.append(lightgbm_score_duration)

fig, axs = plt.subplots(3, sharex=True)

axs[0].plot(n_samples_list, sklearn_scores, label='sklearn')
axs[1].plot(n_samples_list, sklearn_fit_durations, label='sklearn')
axs[2].plot(n_samples_list, sklearn_score_durations, label='sklearn')

if args.lightgbm:
    axs[0].plot(n_samples_list, lightgbm_scores, label='lgbm')
    axs[1].plot(n_samples_list, lightgbm_fit_durations, label='lgbm')
    axs[2].plot(n_samples_list, lightgbm_score_durations, label='lgbm')

for ax in axs:
    ax.set_xscale('log')
    ax.legend(loc='best')
    ax.set_xlabel('n_samples')

axs[0].set_title('scores')
axs[1].set_title('fit duration (s)')
axs[2].set_title('score duration (s)')

title = args.problem
if args.problem == 'classification':
    title += ' n_classes = {}'.format(args.n_classes)
fig.suptitle(title)


plt.tight_layout()
plt.show()
