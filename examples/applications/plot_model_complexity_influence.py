"""
==========================
Model Complexity Influence
==========================

Demonstrate how model complexity influences both prediction accuracy and
computational performance.

The dataset is the boston housing dataset.

"""
from numpy.core.multiarray import count_nonzero
from sklearn.ensemble.gradient_boosting import GradientBoostingRegressor
from sklearn.grid_search import GridSearchCV
from sklearn.linear_model.coordinate_descent import ElasticNet
from sklearn.linear_model.stochastic_gradient import SGDRegressor
from sklearn.pipeline import Pipeline
from sklearn.preprocessing.data import Normalizer

print(__doc__)

# Author: Eustache Diemert <eustache@diemert.fr>
# Forked from the boosting ensemble regression example by
# Peter Prettenhofer <peter.prettenhofer@gmail.com>
# License: BSD 3 clause

import time
import numpy as np
from scipy.sparse import csr_matrix
import pylab as pl

from sklearn import datasets
from sklearn.utils import shuffle
from sklearn.metrics import mean_squared_error
from sklearn.ensemble.forest import ExtraTreesRegressor
from sklearn.svm.classes import NuSVR

from mpl_toolkits.axes_grid1.parasite_axes import host_subplot
from mpl_toolkits.axisartist.axislines import Axes

###############################################################################
# Load data
boston = datasets.load_boston()
X, y = shuffle(boston.data, boston.target, random_state=13)
X = X.astype(np.float32)
offset = int(X.shape[0] * 0.9)
X_train, y_train = X[:offset], y[:offset]
X_test, y_test = X[offset:], y[offset:]
X_test = np.array(X_test)
X_train_normed = Normalizer().fit_transform(X_train)
X_test_normed = Normalizer().fit_transform(X_test)
X_train_sparse = csr_matrix(X_train)
X_test_sparse = csr_matrix(X_test)
X_train_normed_sparse = csr_matrix(X_train_normed)
X_test_normed_sparse = csr_matrix(X_test_normed)

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
        if isinstance(estimator, NuSVR):
            print('using normalized data')
            estimator.fit(X_train_normed, y_train)
        elif isinstance(estimator, ElasticNet):
            print('using sparse data')
            estimator.fit(X_train_normed_sparse, y_train)
        else:
            estimator.fit(X_train, y_train)
        complexity = conf['complexity_computer'](estimator)
        complexities.append(complexity)
        #if isinstance(estimator, ElasticNet):
        #    estimator.sparsify()
        start_time = time.time()
        y_pred = None
        if isinstance(estimator, NuSVR):
            y_pred = estimator.predict(X_test_normed)
        elif isinstance(estimator, ElasticNet):
            y_pred = estimator.predict(X_test_normed_sparse)
        else:
            y_pred = estimator.predict(X_test)
        elapsed_time = time.time() - start_time
        prediction_times.append(elapsed_time)
        mse = mean_squared_error(y_test, y_pred)
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
    y2_label = "Prediction Time (s)"
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
configurations = [{'estimator': ElasticNet,
                   'tuned_params': {'fit_intercept': [True]},
                   'changing_param': 'alpha',
                   'changing_param_values': [1e-6, 1e-5, 1e-4, 1e-3],
                   'complexity_label': 'non_zero coefficients',
                   'complexity_computer':
                       lambda x: count_nonzero(x.coef_)
                   },
                  {'estimator': NuSVR,
                   'tuned_params': {'C': [1e-5, 1e-2, 1, 1e2, 1e3, 1e4],
                                    'gamma': [2**-15, 2**-10, 2**-5,
                                              2**-1, 2, 4, 8]},
                   'changing_param': 'nu',
                   'changing_param_values': [0.1, 0.25, 0.5, 0.75, 0.9],
                   'complexity_label': 'n_support_vectors',
                   'complexity_computer': lambda x: len(x.support_vectors_)
                   },
                  {'estimator': GradientBoostingRegressor,
                   'tuned_params': {'loss': ['ls', 'lad', 'huber',
                                             'quantile']},
                   'changing_param': 'n_estimators',
                   'changing_param_values': [20, 50, 100, 200, 300, 500],
                   'complexity_label': 'n_trees * max_depth',
                   'complexity_computer':
                       lambda x: x.n_estimators*x.max_depth
                   },
                  {'estimator': ExtraTreesRegressor,
                   'tuned_params': {'criterion': ['gini', 'entropy'],
                                    'bootstrap': [True, False],
                                    'n_jobs': 6},
                   'changing_param': 'n_estimators',
                   'changing_param_values': [20, 50, 100, 200, 300, 500],
                   'complexity_label': 'n_trees',
                   'complexity_computer':
                       lambda x: x.n_estimators
                   }]

for conf in configurations:
    clf = conf['estimator']()
    print("Tuning hyperparameters for model: %s" % clf)
    gscv = GridSearchCV(clf, conf['tuned_params'], cv=3, n_jobs=6)
    gscv.fit(X_train, y_train)
    print("Best params: %s" % gscv.best_params_)
    conf['tuned_params'] = gscv.best_params_
    mse_values, prediction_times, complexities = benchmark_influence(conf)
    plot_influence(conf, mse_values, prediction_times, complexities)
