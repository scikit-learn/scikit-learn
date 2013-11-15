"""
==========================
Model Complexity Influence
==========================

Demonstrate how model complexity influences both prediction accuracy and
computational performance.

The dataset is the Boston Housing dataset (regression).

For each class of models we make the model complexity vary through the choice
of relevant model parameters and measure the influence on both computational
performance (latency) and predictive power (MSE).
"""
from sklearn.linear_model.stochastic_gradient import SGDRegressor

print(__doc__)

# Author: Eustache Diemert <eustache@diemert.fr>
# Inspired by the boosting ensemble regression example by
# Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: BSD 3 clause

import time
import numpy as np
from scipy.sparse import csr_matrix
import pylab as pl
from mpl_toolkits.axes_grid1.parasite_axes import host_subplot
from mpl_toolkits.axisartist.axislines import Axes

from sklearn import datasets
from sklearn.utils import shuffle
from sklearn.metrics import mean_squared_error
from sklearn.ensemble.forest import ExtraTreesRegressor
from sklearn.svm.classes import NuSVR
from sklearn.ensemble.gradient_boosting import GradientBoostingRegressor
from sklearn.linear_model.coordinate_descent import ElasticNet
from sklearn import preprocessing
from sklearn.utils.fixes import count_nonzero

###############################################################################
# load data

boston = datasets.load_boston()
np.random.seed(0)


def generate_dataset_versions(normalized=False, sparse=False,
                              enlarge_factor=1, noise=0):
    """Generate sparse/normed/noisy versions of the dataset."""
    X, y = shuffle(boston.data, boston.target)
    #X = X.astype(np.float32)
    #y = y.astype(np.float32)
    if enlarge_factor > 1:
        X = np.repeat(X, repeats=int(enlarge_factor), axis=0)
        y = np.repeat(y, repeats=int(enlarge_factor), axis=0)
    if noise > 0:
        for _ in range(int(noise)):
            y += np.random.rand(*y.shape)
    offset = int(X.shape[0] * 0.8)
    X_train, y_train = X[:offset], y[:offset]
    X_test, y_test = X[offset:], y[offset:]
    X_test = np.array(X_test)
    X_train = np.array(X_train)
    y_test = np.array(y_test)
    y_train = np.array(y_train)
    if normalized:
        min_max_scaler = preprocessing.MinMaxScaler(feature_range=(0, 1))
        min_max_scaler.fit(X_train)
        X_train = min_max_scaler.transform(X_train)
        X_test = min_max_scaler.transform(X_test)
    if sparse:
        X_train = csr_matrix(X_train)
        X_test = csr_matrix(X_test)
    data = {'X_train': X_train, 'X_test': X_test, 'y_train': y_train,
            'y_test': y_test}
    print(data)
    return data


def benchmark_influence(conf):
    """
    Benchmark influence of :changing_param: on both MSE and latency.
    """
    prediction_times = []
    mse_values = []
    complexities = []
    for param_value in conf['changing_param_values']:
        conf['tuned_params'][conf['changing_param']] = param_value
        estimator = conf['estimator'](**conf['tuned_params'])
        print("Benchmarking %s" % estimator)
        estimator.fit(conf['data']['X_train'], conf['data']['y_train'])
        conf['post_fit_hook'](estimator)
        complexity = conf['complexity_computer'](estimator)
        complexities.append(complexity)
        start_time = time.time()
        for _ in range(conf['n_samples']):
            y_pred = estimator.predict(conf['data']['X_test'])
        elapsed_time = (time.time() - start_time) / float(conf['n_samples'])
        prediction_times.append(elapsed_time)
        mse = mean_squared_error(conf['data']['y_test'], y_pred)
        mse_values.append(mse)
        print("Complexity: %d | MSE: %.4f | Pred. Time: %fs\n" % (
            complexity, mse, elapsed_time))
    return mse_values, prediction_times, complexities


def plot_influence(conf, mse_values, prediction_times, complexities):
    """
    Plot influence of model complexity on both accuracy and latency.
    """
    pl.figure(figsize=(12, 6))
    host = host_subplot(111, axes_class=Axes)
    pl.subplots_adjust(right=0.75)
    par1 = host.twinx()
    host.set_xlabel('Model Complexity (%s)' % conf['complexity_label'])
    y1_label = "MSE"
    y2_label = "Time (s)"
    host.set_ylabel(y1_label)
    par1.set_ylabel(y2_label)
    p1, = host.plot(complexities, mse_values, 'b-', label="prediction error")
    p2, = par1.plot(complexities, prediction_times, 'r-',
                    label="latency")
    host.legend(loc='upper right')
    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    pl.title('Influence of Model Complexity - %s' % conf['estimator'].__name__)
    pl.show()

###############################################################################
# main code
configurations = [
                  #{'estimator': NuSVR,
                  # 'tuned_params': {'C': 1e3, 'gamma': 2**-15},
                  # 'changing_param': 'nu',
                  # 'changing_param_values': [0.1, 0.25, 0.5, 0.75, 0.9],
                  # 'complexity_label': 'n_support_vectors',
                  # 'complexity_computer': lambda x: len(x.support_vectors_)
                  # },
                  #{'estimator': GradientBoostingRegressor,
                  # 'tuned_params': {'loss': 'ls'},
                  # 'changing_param': 'n_estimators',
                  # 'changing_param_values': [10, 50, 100, 200, 500],
                  # 'complexity_label': 'n_trees',
                  # 'complexity_computer': lambda x: x.n_estimators
                  # },
                  #{'estimator': ExtraTreesRegressor,
                  # 'tuned_params': {},
                  # 'changing_param': 'n_estimators',
                  # 'changing_param_values': [10, 50, 100, 200, 500],
                  # 'complexity_label': 'n_trees',
                  # 'complexity_computer': lambda x: x.n_estimators
                  # },
                  {'estimator': SGDRegressor,
                   'tuned_params': {'penalty': 'elasticnet',
                                    'alpha': 0.01,
                                    'loss': 'squared_epsilon_insensitive'},
                   'changing_param': 'l1_ratio',
                   'changing_param_values': [0.1, 0.25, 0.5, 0.75, 0.9],
                   'post_fit_hook': lambda x: x.sparsify(),
                   'complexity_label': 'non_zero coefficients',
                   'complexity_computer': lambda x: count_nonzero(x.coef_),
                   'data': generate_dataset_versions(normalized=True,
                                                     sparse=False,
                                                     noise=0),
                   'n_samples': 30
                   },
                  #{'estimator': ElasticNet,
                  # 'tuned_params': {'fit_intercept': True, 'alpha': 0.5},
                  # 'changing_param': 'l1_ratio',
                  # 'changing_param_values': [0.1, 0.25, 0.5, 0.75, 0.9],
                  # 'complexity_label': 'non_zero coefficients',
                  # 'complexity_computer': lambda x: count_nonzero(x.coef_)
                  # }
                ]
for conf in configurations:
    mse_values, prediction_times, complexities = benchmark_influence(conf)
    plot_influence(conf, mse_values, prediction_times, complexities)
