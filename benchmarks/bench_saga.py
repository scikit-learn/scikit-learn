"""Author: Arthur Mensch, Nelle Varoquaux

Benchmarks of sklearn SAGA vs lightning SAGA vs Liblinear. Shows the gain
in using multinomial logistic regression in term of learning time.
"""
import json
import time
import os

from joblib import delayed, Parallel
import matplotlib.pyplot as plt
import numpy as np

from sklearn.datasets import fetch_rcv1, load_iris, load_digits, \
    fetch_20newsgroups_vectorized
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelBinarizer, LabelEncoder
from sklearn.utils.extmath import safe_sparse_dot, softmax


def fit_single(solver, X, y, penalty='l2', single_target=True, C=1,
               max_iter=10, skip_slow=False, dtype=np.float64):
    if skip_slow and solver == 'lightning' and penalty == 'l1':
        print('skip_slowping l1 logistic regression with solver lightning.')
        return

    print('Solving %s logistic regression with penalty %s, solver %s.'
          % ('binary' if single_target else 'multinomial',
             penalty, solver))

    if solver == 'lightning':
        from lightning.classification import SAGAClassifier

    if single_target or solver not in ['sag', 'saga']:
        multi_class = 'ovr'
    else:
        multi_class = 'multinomial'
    X = X.astype(dtype)
    y = y.astype(dtype)
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42,
                                                        stratify=y)
    n_samples = X_train.shape[0]
    n_classes = np.unique(y_train).shape[0]
    test_scores = [1]
    train_scores = [1]
    accuracies = [1 / n_classes]
    times = [0]

    if penalty == 'l2':
        alpha = 1. / (C * n_samples)
        beta = 0
        lightning_penalty = None
    else:
        alpha = 0.
        beta = 1. / (C * n_samples)
        lightning_penalty = 'l1'

    for this_max_iter in range(1, max_iter + 1, 2):
        print('[%s, %s, %s] Max iter: %s' %
              ('binary' if single_target else 'multinomial',
               penalty, solver, this_max_iter))
        if solver == 'lightning':
            lr = SAGAClassifier(loss='log', alpha=alpha, beta=beta,
                                penalty=lightning_penalty,
                                tol=-1, max_iter=this_max_iter)
        else:
            lr = LogisticRegression(solver=solver,
                                    multi_class=multi_class,
                                    C=C,
                                    penalty=penalty,
                                    fit_intercept=False, tol=0,
                                    max_iter=this_max_iter,
                                    random_state=42,
                                    )

        # Makes cpu cache even for all fit calls
        X_train.max()
        t0 = time.clock()

        lr.fit(X_train, y_train)
        train_time = time.clock() - t0

        scores = []
        for (X, y) in [(X_train, y_train), (X_test, y_test)]:
            try:
                y_pred = lr.predict_proba(X)
            except NotImplementedError:
                # Lightning predict_proba is not implemented for n_classes > 2
                y_pred = _predict_proba(lr, X)
            score = log_loss(y, y_pred, normalize=False) / n_samples
            score += (0.5 * alpha * np.sum(lr.coef_ ** 2) +
                      beta * np.sum(np.abs(lr.coef_)))
            scores.append(score)
        train_score, test_score = tuple(scores)

        y_pred = lr.predict(X_test)
        accuracy = np.sum(y_pred == y_test) / y_test.shape[0]
        test_scores.append(test_score)
        train_scores.append(train_score)
        accuracies.append(accuracy)
        times.append(train_time)
    return lr, times, train_scores, test_scores, accuracies


def _predict_proba(lr, X):
    pred = safe_sparse_dot(X, lr.coef_.T)
    if hasattr(lr, "intercept_"):
        pred += lr.intercept_
    return softmax(pred)


def exp(solvers, penalty, single_target,
        n_samples=30000, max_iter=20,
        dataset='rcv1', n_jobs=1, skip_slow=False):
    dtypes_mapping = {
        "float64": np.float64,
        "float32": np.float32,
    }

    if dataset == 'rcv1':
        rcv1 = fetch_rcv1()

        lbin = LabelBinarizer()
        lbin.fit(rcv1.target_names)

        X = rcv1.data
        y = rcv1.target
        y = lbin.inverse_transform(y)
        le = LabelEncoder()
        y = le.fit_transform(y)
        if single_target:
            y_n = y.copy()
            y_n[y > 16] = 1
            y_n[y <= 16] = 0
            y = y_n

    elif dataset == 'digits':
        X, y = load_digits(return_X_y=True)
        if single_target:
            y_n = y.copy()
            y_n[y < 5] = 1
            y_n[y >= 5] = 0
            y = y_n
    elif dataset == 'iris':
        iris = load_iris()
        X, y = iris.data, iris.target
    elif dataset == '20newspaper':
        ng = fetch_20newsgroups_vectorized()
        X = ng.data
        y = ng.target
        if single_target:
            y_n = y.copy()
            y_n[y > 4] = 1
            y_n[y <= 16] = 0
            y = y_n

    X = X[:n_samples]
    y = y[:n_samples]

    out = Parallel(n_jobs=n_jobs, mmap_mode=None)(
        delayed(fit_single)(solver, X, y,
                            penalty=penalty, single_target=single_target,
                            dtype=dtype,
                            C=1, max_iter=max_iter, skip_slow=skip_slow)
        for solver in solvers
        for dtype in dtypes_mapping.values())

    res = []
    idx = 0
    for dtype_name in dtypes_mapping.keys():
        for solver in solvers:
            if not (skip_slow and
                    solver == 'lightning' and
                    penalty == 'l1'):
                lr, times, train_scores, test_scores, accuracies = out[idx]
                this_res = dict(solver=solver, penalty=penalty,
                                dtype=dtype_name,
                                single_target=single_target,
                                times=times, train_scores=train_scores,
                                test_scores=test_scores,
                                accuracies=accuracies)
                res.append(this_res)
            idx += 1

    with open('bench_saga.json', 'w+') as f:
        json.dump(res, f)


def plot(outname=None):
    import pandas as pd
    with open('bench_saga.json', 'r') as f:
        f = json.load(f)
    res = pd.DataFrame(f)
    res.set_index(['single_target'], inplace=True)

    grouped = res.groupby(level=['single_target'])

    colors = {'saga': 'C0', 'liblinear': 'C1', 'lightning': 'C2'}
    linestyles = {"float32": "--", "float64": "-"}
    alpha = {"float64": 0.5, "float32": 1}

    for idx, group in grouped:
        single_target = idx
        fig, axes = plt.subplots(figsize=(12, 4), ncols=4)
        ax = axes[0]

        for scores, times, solver, dtype in zip(group['train_scores'],
                                                group['times'],
                                                group['solver'],
                                                group["dtype"]):
            ax.plot(times, scores, label="%s - %s" % (solver, dtype),
                    color=colors[solver],
                    alpha=alpha[dtype],
                    marker=".",
                    linestyle=linestyles[dtype])
            ax.axvline(times[-1], color=colors[solver],
                       alpha=alpha[dtype],
                       linestyle=linestyles[dtype])
        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Training objective (relative to min)')
        ax.set_yscale('log')

        ax = axes[1]

        for scores, times, solver, dtype in zip(group['test_scores'],
                                                group['times'],
                                                group['solver'],
                                                group["dtype"]):
            ax.plot(times, scores, label=solver, color=colors[solver],
                    linestyle=linestyles[dtype],
                    marker=".",
                    alpha=alpha[dtype])
            ax.axvline(times[-1], color=colors[solver],
                       alpha=alpha[dtype],
                       linestyle=linestyles[dtype])

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Test objective (relative to min)')
        ax.set_yscale('log')

        ax = axes[2]
        for accuracy, times, solver, dtype in zip(group['accuracies'],
                                                  group['times'],
                                                  group['solver'],
                                                  group["dtype"]):
            ax.plot(times, accuracy, label="%s - %s" % (solver, dtype),
                    alpha=alpha[dtype],
                    marker=".",
                    color=colors[solver], linestyle=linestyles[dtype])
            ax.axvline(times[-1], color=colors[solver],
                       alpha=alpha[dtype],
                       linestyle=linestyles[dtype])

        ax.set_xlabel('Time (s)')
        ax.set_ylabel('Test accuracy')
        ax.legend()
        name = 'single_target' if single_target else 'multi_target'
        name += '_%s' % penalty
        plt.suptitle(name)
        if outname is None:
            outname = name + '.png'
        fig.tight_layout()
        fig.subplots_adjust(top=0.9)

        ax = axes[3]
        for scores, times, solver, dtype in zip(group['train_scores'],
                                                group['times'],
                                                group['solver'],
                                                group["dtype"]):
            ax.plot(np.arange(len(scores)),
                    scores, label="%s - %s" % (solver, dtype),
                    marker=".",
                    alpha=alpha[dtype],
                    color=colors[solver], linestyle=linestyles[dtype])

        ax.set_yscale("log")
        ax.set_xlabel('# iterations')
        ax.set_ylabel('Objective function')
        ax.legend()

        plt.savefig(outname)


if __name__ == '__main__':
    solvers = ['saga', 'liblinear', 'lightning']
    penalties = ['l1', 'l2']
    n_samples = [100000, 300000, 500000, 800000, None]
    single_target = True
    for penalty in penalties:
        for n_sample in n_samples:
            exp(solvers, penalty, single_target,
                n_samples=n_sample, n_jobs=1,
                dataset='rcv1', max_iter=10)
            if n_sample is not None:
                outname = "figures/saga_%s_%d.png" % (penalty, n_sample)
            else:
                outname = "figures/saga_%s_all.png" % (penalty,)
            try:
                os.makedirs("figures")
            except OSError:
                pass
            plot(outname)
