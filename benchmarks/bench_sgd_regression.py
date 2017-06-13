"""
Benchmark for SGD regression

Compares SGD regression against coordinate descent and Ridge
on synthetic data.
"""

print(__doc__)

# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt

import gc

from time import time

from sklearn.linear_model import Ridge, SGDRegressor, ElasticNet
from sklearn.metrics import mean_squared_error
from sklearn.datasets.samples_generator import make_regression

if __name__ == "__main__":
    list_n_samples = np.linspace(100, 10000, 5).astype(np.int)
    list_n_features = [10, 100, 1000]
    n_test = 1000
    noise = 0.1
    alpha = 0.01
    sgd_results = np.zeros((len(list_n_samples), len(list_n_features), 2))
    elnet_results = np.zeros((len(list_n_samples), len(list_n_features), 2))
    ridge_results = np.zeros((len(list_n_samples), len(list_n_features), 2))
    asgd_results = np.zeros((len(list_n_samples), len(list_n_features), 2))
    for i, n_train in enumerate(list_n_samples):
        for j, n_features in enumerate(list_n_features):
            X, y, coef = make_regression(
                n_samples=n_train + n_test, n_features=n_features,
                noise=noise, coef=True)

            X_train = X[:n_train]
            y_train = y[:n_train]
            X_test = X[n_train:]
            y_test = y[n_train:]

            print("=======================")
            print("Round %d %d" % (i, j))
            print("n_features:", n_features)
            print("n_samples:", n_train)

            # Shuffle data
            idx = np.arange(n_train)
            np.random.seed(13)
            np.random.shuffle(idx)
            X_train = X_train[idx]
            y_train = y_train[idx]

            std = X_train.std(axis=0)
            mean = X_train.mean(axis=0)
            X_train = (X_train - mean) / std
            X_test = (X_test - mean) / std

            std = y_train.std(axis=0)
            mean = y_train.mean(axis=0)
            y_train = (y_train - mean) / std
            y_test = (y_test - mean) / std

            gc.collect()
            print("- benchmarking ElasticNet")
            clf = ElasticNet(alpha=alpha, l1_ratio=0.5, fit_intercept=False)
            tstart = time()
            clf.fit(X_train, y_train)
            elnet_results[i, j, 0] = mean_squared_error(clf.predict(X_test),
                                                       y_test)
            elnet_results[i, j, 1] = time() - tstart

            gc.collect()
            print("- benchmarking SGD")
            n_iter = np.ceil(10 ** 4.0 / n_train)
            clf = SGDRegressor(alpha=alpha / n_train, fit_intercept=False,
                               n_iter=n_iter, learning_rate="invscaling",
                               eta0=.01, power_t=0.25)

            tstart = time()
            clf.fit(X_train, y_train)
            sgd_results[i, j, 0] = mean_squared_error(clf.predict(X_test),
                                                     y_test)
            sgd_results[i, j, 1] = time() - tstart

            gc.collect()
            print("n_iter", n_iter)
            print("- benchmarking A-SGD")
            n_iter = np.ceil(10 ** 4.0 / n_train)
            clf = SGDRegressor(alpha=alpha / n_train, fit_intercept=False,
                               n_iter=n_iter, learning_rate="invscaling",
                               eta0=.002, power_t=0.05,
                               average=(n_iter * n_train // 2))

            tstart = time()
            clf.fit(X_train, y_train)
            asgd_results[i, j, 0] = mean_squared_error(clf.predict(X_test),
                                                       y_test)
            asgd_results[i, j, 1] = time() - tstart

            gc.collect()
            print("- benchmarking RidgeRegression")
            clf = Ridge(alpha=alpha, fit_intercept=False)
            tstart = time()
            clf.fit(X_train, y_train)
            ridge_results[i, j, 0] = mean_squared_error(clf.predict(X_test),
                                                       y_test)
            ridge_results[i, j, 1] = time() - tstart

    # Plot results
    i = 0
    m = len(list_n_features)
    plt.figure('scikit-learn SGD regression benchmark results',
               figsize=(5 * 2, 4 * m))
    for j in range(m):
        plt.subplot(m, 2, i + 1)
        plt.plot(list_n_samples, np.sqrt(elnet_results[:, j, 0]),
                 label="ElasticNet")
        plt.plot(list_n_samples, np.sqrt(sgd_results[:, j, 0]),
                 label="SGDRegressor")
        plt.plot(list_n_samples, np.sqrt(asgd_results[:, j, 0]),
                 label="A-SGDRegressor")
        plt.plot(list_n_samples, np.sqrt(ridge_results[:, j, 0]),
                 label="Ridge")
        plt.legend(prop={"size": 10})
        plt.xlabel("n_train")
        plt.ylabel("RMSE")
        plt.title("Test error - %d features" % list_n_features[j])
        i += 1

        plt.subplot(m, 2, i + 1)
        plt.plot(list_n_samples, np.sqrt(elnet_results[:, j, 1]),
                 label="ElasticNet")
        plt.plot(list_n_samples, np.sqrt(sgd_results[:, j, 1]),
                 label="SGDRegressor")
        plt.plot(list_n_samples, np.sqrt(asgd_results[:, j, 1]),
                 label="A-SGDRegressor")
        plt.plot(list_n_samples, np.sqrt(ridge_results[:, j, 1]),
                 label="Ridge")
        plt.legend(prop={"size": 10})
        plt.xlabel("n_train")
        plt.ylabel("Time [sec]")
        plt.title("Training time - %d features" % list_n_features[j])
        i += 1

    plt.subplots_adjust(hspace=.30)

    plt.show()
