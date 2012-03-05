"""
Benchmark script to bench scikit-learn GradientBoosting vs. R's gbm package .

NOTE::

make sure you run
$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib64/R/lib

"""

import numpy as np
import rpy2

from time import time
from functools import wraps

from sklearn import datasets
from sklearn.utils import shuffle
from sklearn.utils import check_random_state
from sklearn.externals import joblib

from rpy2.robjects.numpy2ri import numpy2ri
from rpy2.robjects.packages import importr

import pylab as pl

gbm = importr('gbm')


def repeat(n_repetitions=10):
    def wrap(f):
        def wrapper(*args, **kargs):
            scores = []
            for i in range(n_repetitions):
                scores.append(f(*args, random_state=i, **kargs))
            scores = np.array(scores)
            return scores.mean(axis=0), scores.std(axis=0)
        return wraps(f)(wrapper)
    return wrap


# ignore overflows due to exp
#np.seterr(invalid='print', under='print', divide='print', over='ignore')

classification_params = {"distribution": "bernoulli", "shrinkage": 0.6,
                         "n.tree": 250, "bag.fraction": 1.0, "verbose": False,
                         "n.minobsinnode": 4, "interaction.depth": 1, "offset": False}

@repeat()
def bench_random_gaussian(random_state=None):
    rs = check_random_state(random_state)
    shape = (12000, 10)
    X = rs.normal(size=shape).reshape(shape)
    y = ((X ** 2.0).sum(axis=1) > 9.34).astype(np.float64)

    X_train, X_test = X[:2000], X[2000:]
    y_train, y_test = y[:2000], y[2000:]

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **classification_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":classification_params["n.tree"]})
    test_time = time() - t0
    pred = (np.array(pred) >= 0.0).astype(np.float64)
    error_rate = np.mean(pred != y_test)
    #trees = [gbm.pretty_gbm_tree(model, **{'i.tree': i+1}) for i in range(classification_params['n.tree'])]
    #avg_splits = np.mean([sum(1 for n in t[2] if n != -1) for t in trees])
    #print "Avg. splits: %.4f" % avg_splits
    return error_rate, train_time, test_time


@repeat()
def bench_spam(random_state=None):
    X = np.loadtxt("/home/pprett/corpora/spam/spambase.data", delimiter=",")
    y = X[:, -1].ravel()
    X = X[:, :-1]
    f = open("/home/pprett/corpora/spam/spambase.names")
    feature_names = np.array([l.split(":")[0] for l in f])

    X, y = shuffle(X, y, random_state=random_state)
    X_test, y_test = X[:1536], y[:1536]
    X_train, y_train = X[1536:], y[1536:]

    y_train[y_train == -1.0] = 0
    y_test[y_test == -1.0] = 0

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **classification_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":classification_params["n.tree"]})
    test_time = time() - t0
    pred = (np.array(pred) >= 0.0).astype(np.float64)
    error_rate = np.mean(pred != y_test)

    return error_rate, train_time, test_time


@repeat()
def bench_madelon(random_state=None):
    X_train = np.loadtxt("/home/pprett/corpora/madelon/madelon_train.data")
    y_train = np.loadtxt("/home/pprett/corpora/madelon/madelon_train.labels")
    X_test = np.loadtxt("/home/pprett/corpora/madelon/madelon_valid.data")
    y_test = np.loadtxt("/home/pprett/corpora/madelon/madelon_valid.labels")

    y_train[y_train == -1] = 0
    y_test[y_test == -1] = 0

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **classification_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":classification_params["n.tree"]})
    test_time = time() - t0
    pred = (np.array(pred) >= 0.0).astype(np.float64)
    error_rate = np.mean(pred != y_test)

    return error_rate, train_time, test_time


@repeat()
def bench_arcene(random_state=None):
    X_train = np.loadtxt("/home/pprett/corpora/arcene/arcene_train.data")
    y_train = np.loadtxt("/home/pprett/corpora/arcene/arcene_train.labels")
    X_test = np.loadtxt("/home/pprett/corpora/arcene/arcene_valid.data")
    y_test = np.loadtxt("/home/pprett/corpora/arcene/arcene_valid.labels")

    y_train[y_train == -1.0] = 0
    y_test[y_test == -1.0] = 0

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **classification_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":classification_params["n.tree"]})
    test_time = time() - t0
    pred = (np.array(pred) >= 0.0).astype(np.float64)
    error_rate = np.mean(pred != y_test)

    return error_rate, train_time, test_time


regression_params = {"distribution": "gaussian", "shrinkage": 0.1,
                     "n.tree": 250, "bag.fraction": 1.0, "verbose": False,
                     "n.minobsinnode": 1, "interaction.depth": 1,
                     "offset": False}

@repeat()
def bench_boston(random_state=None):
    boston = datasets.load_boston()
    X, y = shuffle(boston.data, boston.target, random_state=random_state)
    X = X.astype(np.float32)
    y = y.astype(np.float32)
    offset = int(X.shape[0] * 0.9)
    X_train = X[:offset]
    y_train = y[:offset]
    X_test = X[offset:]
    y_test = y[offset:]

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **regression_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":regression_params["n.tree"]})
    test_time = time() - t0
    pred = np.array(pred, dtype=np.float64)

    mse = np.mean((pred - y_test) ** 2.0)
    return mse, train_time, test_time


@repeat()
def bench_friedman1(random_state=None):
    X, y = datasets.make_friedman1(n_samples=1200,
                                   random_state=random_state, noise=1.0)
    X = X.astype(np.float32)
    y = y.astype(np.float32)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **regression_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":regression_params["n.tree"]})
    test_time = time() - t0
    pred = np.array(pred, dtype=np.float64)
    mse = np.mean((pred - y_test) ** 2.0)
    return mse, train_time, test_time


@repeat()
def bench_friedman2(random_state=None):
    X, y = datasets.make_friedman2(n_samples=1200, random_state=random_state)
    X = X.astype(np.float32)
    y = y.astype(np.float32)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **regression_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":regression_params["n.tree"]})
    test_time = time() - t0
    pred = np.array(pred, dtype=np.float64)
    mse = np.mean((pred - y_test) ** 2.0)
    return mse, train_time, test_time


@repeat()
def bench_friedman3(random_state=None):
    X, y = datasets.make_friedman3(n_samples=1200, random_state=random_state)
    X = X.astype(np.float32)
    y = y.astype(np.float32)
    X_train, y_train = X[:200], y[:200]
    X_test, y_test = X[200:], y[200:]

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **regression_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":regression_params["n.tree"]})
    test_time = time() - t0
    pred = np.array(pred, dtype=np.float64)
    mse = np.mean((pred - y_test) ** 2.0)
    return mse, train_time, test_time


@repeat(1)
def bench_yahoo_ltrc(random_state=None):
    data = joblib.load('/home/pprett/corpora/yahoo-ltrc-2010/data/set1.pkl')

    X_train = data['X_train']
    y_train = data['y_train']

    X_test = data['X_test']
    y_test = data['y_test']

    X_train = numpy2ri(X_train)
    X_test = numpy2ri(X_test)
    y_train = numpy2ri(y_train)

    t0 = time()
    model = gbm.gbm_fit(X_train, y_train, **regression_params)
    train_time = time() - t0
    t0 = time()
    pred = gbm.predict_gbm(model, X_test,
                           **{"n.tree":regression_params["n.tree"]})
    test_time = time() - t0
    pred = np.array(pred, dtype=np.float64)
    mse = np.mean((pred - y_test) ** 2.0)
    return mse, train_time, test_time


gbm_results = {
    "Example 10.2": bench_random_gaussian(),
    "Spam": bench_spam(),
    "Madelon": bench_madelon(),
    "Arcene": bench_arcene(),
    "Boston": bench_boston(),
    "Friedman#1": bench_friedman1(),
    "Friedman#2": bench_friedman2(),
    "Friedman#3": bench_friedman3(),
    "YahooLTRC": bench_yahoo_ltrc(),
    }

from pprint import pprint
pprint(gbm_results)


def plot_against(gbrt_results, gbm_results, clf=True):

    if clf:
        datasets = [
            "Example 10.2",
            "Spam",
            "Madelon",
            "Arcene",
            ]
    else:
        datasets = [
            "Boston",
            "Friedman#1",
            "Friedman#2",
            "Friedman#3",
        ]

    width = 0.35
    ind = np.arange(len(datasets))

    fig = pl.figure()

    for i, measure in enumerate(['error rate', 'train time', 'test time']):
        ax = fig.add_subplot(1, 3, i + 1)

        gbrt_means = np.array([gbrt_results[ds][0][i] for ds in datasets])
        gbrt_std = np.array([gbrt_results[ds][1][i] for ds in datasets])

        gbm_means = np.array([gbm_results[ds][0][i] for ds in datasets])
        gbm_std = np.array([gbm_results[ds][1][i] for ds in datasets])

        if not clf:
            gbrt_means *= (1.0 / gbm_means)
            gbrt_var = gbrt_std ** 2.0
            gbrt_var *= (1.0 / gbm_means) ** 2.0
            gbrt_std = np.sqrt(gbrt_var)

            gbm_var = gbm_std ** 2.0
            gbm_var *= (1.0 / gbm_means) ** 2.0
            gbm_std = np.sqrt(gbm_var)
            gbm_means *= (1.0 / gbm_means)

            rects1 = ax.bar(ind, gbrt_means, width, color='r', yerr=gbrt_std)
            rects2 = ax.bar(ind + width, gbm_means, width, color='y', yerr=gbm_std)

        else:
            rects1 = ax.bar(ind, gbrt_means, width, color='r', yerr=gbrt_std)
            rects2 = ax.bar(ind + width, gbm_means, width, color='y', yerr=gbm_std)

        ax.set_ylabel(measure)
        ax.set_title("%s: Sklearn vs. R" % measure.title())
        ax.set_xticks(ind+width)
        ax.set_xticklabels(datasets, rotation=45, ha='right')

        ax.legend((rects1[0], rects2[0]), ('Sklearn', 'R (gbm)'), loc='upper left')
