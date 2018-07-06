"""
Benchmarks of sklearn solver in LogisticRegression.
"""

# Author: Tom Dupre la Tour
import time
from os.path import expanduser

import matplotlib.pyplot as plt
import scipy.sparse as sp  # noqa
import numpy as np
import pandas as pd

from sklearn.datasets import fetch_mldata
from sklearn.datasets import fetch_rcv1, load_iris, load_digits
from sklearn.datasets import fetch_20newsgroups_vectorized
from sklearn.exceptions import ConvergenceWarning
from sklearn.externals.joblib import delayed, Parallel, Memory
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model.logistic import _multinomial_loss
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import MinMaxScaler  # noqa
from sklearn.utils.testing import ignore_warnings
from sklearn.utils import shuffle


def get_loss(coefs, intercepts, X, y, C, multi_class, penalty):
    if multi_class == 'ovr':
        if np.array(intercepts).ndim == 0 and intercepts == 0:
            intercepts = np.zeros(coefs.shape[0])
        loss = 0
        for ii, (coef, intercept) in enumerate(zip(coefs, intercepts)):
            y_bin = y.copy()
            y_bin[y == ii] = 1
            y_bin[y != ii] = -1
            loss += np.sum(
                np.log(1. + np.exp(-y_bin * (X.dot(coef) + intercept))))

            if penalty == 'l2':
                loss += 0.5 * coef.dot(coef) / C
            else:
                loss += np.sum(np.abs(coef)) / C
    else:
        coefs_and_intercept = np.vstack((coefs.T, intercepts.T)).T.ravel()
        lbin = LabelBinarizer()
        Y_multi = lbin.fit_transform(y)
        if Y_multi.shape[1] == 1:
            Y_multi = np.hstack([1 - Y_multi, Y_multi])
        loss, _, _ = _multinomial_loss(coefs_and_intercept, X, Y_multi, 0,
                                       np.ones(X.shape[0]))
        coefs = coefs.ravel()
        if penalty == 'l2':
            loss += 0.5 * coefs.dot(coefs) / C
        else:
            loss += np.sum(np.abs(coefs)) / C

    loss /= X.shape[0]

    return loss


def fit_single(solver, X, y, X_shape, dataset, penalty='l2',
               multi_class='multinomial', C=1, max_iter=10):
    assert X.shape == X_shape

    # if not sp.issparse(X):
    #     X = MinMaxScaler().fit_transform(X)

    X_train, X_test, y_train, y_test = train_test_split(
        X, y, random_state=42, stratify=y)
    train_scores, train_losses, test_scores, times = [], [], [], []

    if solver == 'newton-cg':
        max_iter /= 2

    n_repeats = None
    max_iter_range = np.unique(np.int_(np.logspace(0, np.log10(max_iter), 10)))
    for this_max_iter in max_iter_range:
        msg = ('[%s, %s, %s, %s] Max iter: %s' % (dataset, multi_class, solver,
                                                  penalty, this_max_iter))
        lr = LogisticRegression(solver=solver, multi_class=multi_class, C=C,
                                penalty=penalty, fit_intercept=False,
                                tol=1e-24, max_iter=this_max_iter,
                                random_state=42, intercept_scaling=10000)
        t0 = time.clock()
        try:
            if penalty == 'l1' and multi_class == 'multinomial':
                raise ValueError('skip as only saga is available.')

            with ignore_warnings(category=ConvergenceWarning):
                # first time for timing
                if n_repeats is None:
                    t0 = time.clock()
                    lr.fit(X_train, y_train)
                    max_iter_duration = max_iter * (time.clock() - t0)
                    n_repeats = max(1, int(1. / max_iter_duration))

                t0 = time.clock()
                for _ in range(n_repeats):
                    lr.fit(X_train, y_train)
                train_time = (time.clock() - t0) / n_repeats
                print('%s (repeat=%d)' % (msg, n_repeats))

        except ValueError:
            train_score = np.nan
            train_loss = np.nan
            test_score = np.nan
            train_time = np.nan
            print('%s (skipped)' % (msg, ))
            continue

        train_loss = get_loss(lr.coef_, lr.intercept_, X_train, y_train, C,
                              multi_class, penalty)
        train_score = lr.score(X_train, y_train)
        test_score = lr.score(X_test, y_test)

        train_scores.append(train_score)
        train_losses.append(train_loss)
        test_scores.append(test_score)
        times.append(train_time)

    return (solver, penalty, dataset, multi_class, times, train_losses,
            train_scores, test_scores)


def load_dataset(dataset, n_samples_max):
    if dataset == 'rcv1':
        rcv1 = fetch_rcv1()
        X = rcv1.data
        y = rcv1.target

        # take only 3 categories (CCAT, ECAT, MCAT)
        y = y[:, [1, 4, 10]].astype(np.float64)
        # remove samples that have more than one category
        mask = np.asarray(y.sum(axis=1) == 1).ravel()
        y = y[mask, :].indices
        X = X[mask, :]

    elif dataset == 'mnist':
        mnist = fetch_mldata('MNIST original')
        X, y = shuffle(mnist.data, mnist.target, random_state=42)
        X = X.astype(np.float64)

    elif dataset == 'digits':
        digits = load_digits()
        X, y = digits.data, digits.target

    elif dataset == 'iris':
        iris = load_iris()
        X, y = iris.data, iris.target

    elif dataset == '20news':
        ng = fetch_20newsgroups_vectorized()
        X = ng.data
        y = ng.target

    X = X[:n_samples_max]
    y = y[:n_samples_max]

    return X, y


def run(solvers, penalties, multi_classes, n_samples_max, max_iter, datasets,
        n_jobs):
    mem = Memory(cachedir=expanduser('~/cache'), verbose=0)

    results = []
    for dataset in datasets:
        for multi_class in multi_classes:
            X, y = load_dataset(dataset, n_samples_max)

            cached_fit = mem.cache(fit_single, ignore=['X'])
            cached_fit = fit_single

            out = Parallel(n_jobs=n_jobs, mmap_mode=None)(delayed(cached_fit)(
                solver, X, y, X.shape, dataset=dataset, penalty=penalty,
                multi_class=multi_class, C=1, max_iter=max_iter)
                for solver in solvers
                for penalty in penalties)  # yapf: disable

            results.extend(out)

            columns = ("solver penalty dataset multi_class times "
                       "train_losses train_scores test_scores").split()
            results_df = pd.DataFrame(out, columns=columns)
            plot(results_df)


def plot(res):
    res.set_index(['dataset', 'multi_class', 'penalty'], inplace=True)

    grouped = res.groupby(level=['dataset', 'multi_class', 'penalty'])

    colors = {
        'sag': 'red',
        'saga': 'orange',
        'liblinear': 'blue',
        'lbfgs': 'green',
        'newton-cg': 'darkviolet',
        'auto': 'black',
    }

    for idx, group in grouped:
        dataset, multi_class, penalty = idx
        fig = plt.figure(figsize=(12, 4))

        # -----------------------
        ax = fig.add_subplot(131)
        train_losses = group['train_losses']
        tmp = np.sort(np.concatenate(train_losses.values))
        if tmp.size == 0:
            plt.close(fig)
            continue
        ref = 2 * tmp[0] - tmp[1]

        for losses, times, solver in zip(group['train_losses'], group['times'],
                                         group['solver']):
            losses = losses - ref
            linestyle = ':' if solver == 'auto' else '-'
            ax.plot(times, losses, label=solver, color=colors[solver],
                    linestyle=linestyle, marker='.')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Training objective (relative to min)')
        ax.set_yscale('log')

        # -----------------------
        ax = fig.add_subplot(132)

        for train_score, times, solver in zip(group['train_scores'],
                                              group['times'], group['solver']):
            linestyle = ':' if solver == 'auto' else '-'
            ax.plot(times, train_score, label=solver, color=colors[solver],
                    linestyle=linestyle, marker='.')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Train score')

        # -----------------------
        ax = fig.add_subplot(133)

        for test_score, times, solver in zip(group['test_scores'],
                                             group['times'], group['solver']):
            linestyle = ':' if solver == 'auto' else '-'
            ax.plot(times, test_score, label=solver, color=colors[solver],
                    linestyle=linestyle, marker='.')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Test score')
        ax.legend()

        # -----------------------
        name = '%s_%s_%s' % (multi_class, penalty, dataset)
        plt.suptitle(name)
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)
        plt.savefig('figures/' + name + '.png')
        plt.close(fig)
        print('SAVED: ' + name)


if __name__ == '__main__':
    n_jobs = 3
    max_iter = 50
    solvers = ['liblinear', 'saga', 'sag', 'lbfgs', 'newton-cg', 'auto']
    penalties = ['l2', 'l1']
    multi_classes = ['multinomial', 'ovr']
    datasets = ['iris', 'digits', '20news', 'rcv1', 'mnist']

    run(solvers, penalties, multi_classes, n_samples_max=None, n_jobs=n_jobs,
        datasets=datasets, max_iter=max_iter)
