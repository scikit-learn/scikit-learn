"""
This benchmark compares the predictions of a RandomForestRegressor
when the parallelization is disabled. The benchmark depends on
four parameters:

* number of estimators (-e)
* number of observations in the batch (-n)
* number of features (-f)
* number of repetitions (-r)
* use assume_finite (-a)

The first execution fails after saving model to make
sure the training is not part of the benchmark.
"""
import argparse
import time
import pickle
import os
import numpy as np
from sklearn import config_context
from sklearn.datasets import make_regression
from sklearn.ensemble import RandomForestRegressor


def build_model(e, n, f):
    filename = "RandomForestRegressor-e%d-n%d-f%d.pkl" % (e, n, f)
    if os.path.exists(filename):
        with open(filename, "rb") as f:
            return pickle.load(f)

    nt = 10000
    X, y = make_regression(n_samples=nt + n, n_features=f, n_informative= f // 2,
                           n_targets=1, random_state=1)
    X_train, X_test = X[:nt], X[nt:]
    y_train, y_test = y[:nt], y[nt:]
    rf = RandomForestRegressor(n_estimators=e, random_state=1)
    rf.fit(X, y)
    data = dict(model=rf, data=X_test.astype(np.float32))
    with open(filename, "wb") as f:
        pickle.dump(data, f)
    raise RuntimeError("Data and data are cached. Run the script "
                       "again to measure the performance.")


def _run_predict(model, X, repeat, parallel_predict, check_input):
    for r in range(repeat):
        model.predict(X, parallel_predict=parallel_predict,
                      check_input=check_input)


def pp_true(model, X, repeat):
    _run_predict(model, X, repeat, True, True)


def pp_false(model, X, repeat):
    _run_predict(model, X, repeat, False, True)


def pp_false_nocheck(model, X, repeat):
    _run_predict(model, X, repeat, False, False)


def benchmark(model, X, repeat):
    begin = time.perf_counter()
    pp_true(model, X, repeat)
    end1 = time.perf_counter()
    pp_false(model, X, repeat)
    end2 = time.perf_counter()
    pp_false_nocheck(model, X, repeat)
    end3 = time.perf_counter()
    r = repeat
    e = len(model.estimators_)
    n = X.shape[0]
    f = X.shape[1]
    print("parallel=T check=T: r=%d e=%d n=%d f=%d time=%f" % (r, e, n, f, end1 - begin))
    print("parallel=F check=T: r=%d e=%d n=%d f=%d time=%f" % (r, e, n, f, end2 - end1))
    print("parallel=F check=F: r=%d e=%d n=%d f=%d time=%f" % (r, e, n, f, end3 - end2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', type=int, default=100)
    parser.add_argument('-n', type=int, default=1)
    parser.add_argument('-f', type=int, default=10)
    parser.add_argument('-r', type=int, default=1000)
    parser.add_argument('-a', type=int, default=1)
    args = parser.parse_args()

    model_data = build_model(args.e, args.n, args.f)

    if args.a:
        with config_context(assume_finite=True):
            benchmark(model_data['model'], model_data['data'], args.r)
    else:
        benchmark(model_data['model'], model_data['data'], args.r)


main()

# py-spy record --native --function --rate=10 -o profile.svg -- python bench_random_forest_parallel_predict.py -e 100 -n 1 -f 10 -r 1000
