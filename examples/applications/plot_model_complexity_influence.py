"""
==========================
Model Complexity Influence
==========================

Demonstrate how model complexity influences both prediction accuracy and
computational performance.

We will be using two datasets:
    - :ref:`diabetes_dataset` for regression.
      This dataset consists of 10 measurements taken from diabetes patients.
      The task is to predict disease
      progression.
    - :ref:`20newsgroups_dataset` for classification. This dataset consists of
      newsgroup posts. The task is to predict on which topic (out of 20 topics)
      the post is written about.

We will model the complexity influence on three different estimators:
    - :class:`~sklearn.linear_model.SGDClassifier` (for classification data)
      which implements stochastic gradient descent learning

    - :class:`~sklearn.svm.NuSVR` (for regression data) which implements
      Nu support vector regression

    - :class:`~sklearn.ensemble.GradientBoostingRegressor` (for regression
      data) which builds an additive
      model in a forward stage-wise fashion


We make the model complexity vary through the choice of relevant model
parameters in each of our selected models. Next, we will measure the influence
on both computational performance (latency) and predictive power (MSE or
Hamming Loss).

"""

print(__doc__)

# Authors: Eustache Diemert <eustache@diemert.fr>
#          Maria Telenczuk <https://github.com/maikia>
# License: BSD 3 clause

import time
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.parasite_axes import host_subplot
from mpl_toolkits.axisartist.axislines import Axes

from sklearn import datasets
from sklearn.utils import shuffle
from sklearn.metrics import mean_squared_error
from sklearn.svm import NuSVR
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.linear_model import SGDClassifier
from sklearn.metrics import hamming_loss


# Initialize random generator
np.random.seed(0)

##############################################################################
# Load the data
# -------------
#
# First we load both datasets.
#
# Note 1: We are using :func:`~sklearn.datasets.fetch_20newsgroups_vectorized`
#  to download 20
# newsgroups dataset. It returns ready-to-use features.
#
# Note 2: ``X`` of the 20 newsgropus dataset is a sparse matrix while ``X`` of
# diabetes dataset is a numpy array.
#


def generate_data(case):
    """Generate regression/classification data."""
    if case == 'regression':
        X, y = datasets.load_diabetes(return_X_y=True)
    elif case == 'classification':
        X, y = datasets.fetch_20newsgroups_vectorized(subset='all',
                                                      return_X_y=True)
    X, y = shuffle(X, y)
    offset = int(X.shape[0] * 0.8)
    X_train, y_train = X[:offset], y[:offset]
    X_test, y_test = X[offset:], y[offset:]

    data = {'X_train': X_train, 'X_test': X_test, 'y_train': y_train,
            'y_test': y_test}
    return data


regression_data = generate_data('regression')
classification_data = generate_data('classification')


##############################################################################
# Benchmark influence
# -------------------
# Next, we can calculate the influence of the parameters on the given
# estimator. In each lap we will set the estimator with the new value of
# ``changing_param`` and we will be collecting the prediction times, prediction
# performance and complexities to see how those changes affect the estimator.
# We will calculate the complexity using ``complexity_computer`` passed as a
# parameter.
#


def benchmark_influence(conf):
    """
    Benchmark influence of :changing_param: on both MSE and latency.
    """
    prediction_times = []
    prediction_powers = []
    complexities = []
    for param_value in conf['changing_param_values']:
        conf['tuned_params'][conf['changing_param']] = param_value
        estimator = conf['estimator'](**conf['tuned_params'])

        print("Benchmarking %s" % estimator)
        estimator.fit(conf['data']['X_train'], conf['data']['y_train'])
        conf['postfit_hook'](estimator)
        complexity = conf['complexity_computer'](estimator)
        complexities.append(complexity)
        start_time = time.time()
        for _ in range(conf['n_samples']):
            y_pred = estimator.predict(conf['data']['X_test'])
        elapsed_time = (time.time() - start_time) / float(conf['n_samples'])
        prediction_times.append(elapsed_time)
        pred_score = conf['prediction_performance_computer'](
            conf['data']['y_test'], y_pred)
        prediction_powers.append(pred_score)
        print("Complexity: %d | %s: %.4f | Pred. Time: %fs\n" % (
            complexity, conf['prediction_performance_label'], pred_score,
            elapsed_time))
    return prediction_powers, prediction_times, complexities


##############################################################################
# Choose parameters
# -----------------------------
#
# We choose the parameters for each of our estimators by making
# a dictionary with all the necessary values.
# ``changing_param`` is the name of the parameter which will vary in each
# estimator.
# Complexity will be defined by the ``complexity_label`` and calculated using
# 'complexity_computer'
# Also note that depending on the estimator type we are passing
# different data.
#

def _count_nonzero_coefficients(estimator):
    a = estimator.coef_.toarray()
    return np.count_nonzero(a)


configurations = [
    {'estimator': SGDClassifier,
     'tuned_params': {'penalty': 'elasticnet', 'alpha': 0.001, 'loss':
                      'modified_huber', 'fit_intercept': True, 'tol': 1e-3},
     'changing_param': 'l1_ratio',
     'changing_param_values': [0.25, 0.5, 0.75, 0.9],
     'complexity_label': 'non_zero coefficients',
     'complexity_computer': _count_nonzero_coefficients,
     'prediction_performance_computer': hamming_loss,
     'prediction_performance_label': 'Hamming Loss (Misclassification Ratio)',
     'postfit_hook': lambda x: x.sparsify(),
     'data': classification_data,
     'n_samples': 30},
    {'estimator': NuSVR,
     'tuned_params': {'C': 1e3, 'gamma': 2 ** -15},
     'changing_param': 'nu',
     'changing_param_values': [0.1, 0.25, 0.5, 0.75, 0.9],
     'complexity_label': 'n_support_vectors',
     'complexity_computer': lambda x: len(x.support_vectors_),
     'data': regression_data,
     'postfit_hook': lambda x: x,
     'prediction_performance_computer': mean_squared_error,
     'prediction_performance_label': 'MSE',
     'n_samples': 30},
    {'estimator': GradientBoostingRegressor,
     'tuned_params': {'loss': 'ls'},
     'changing_param': 'n_estimators',
     'changing_param_values': [10, 50, 100, 200, 500],
     'complexity_label': 'n_trees',
     'complexity_computer': lambda x: x.n_estimators,
     'data': regression_data,
     'postfit_hook': lambda x: x,
     'prediction_performance_computer': mean_squared_error,
     'prediction_performance_label': 'MSE',
     'n_samples': 30},
]


##############################################################################
# Run the code and plot the results
# ---------------------------------
# We are ready to run all our functions by looping though the configurations we
# set previously and then plotting the influence of varying parameters on the
# complexity and the latency of each model.
#


def plot_influence(conf, mse_values, prediction_times, complexities):
    """
    Plot influence of model complexity on both accuracy and latency.
    """
    plt.figure(figsize=(12, 6))
    host = host_subplot(111, axes_class=Axes)
    plt.subplots_adjust(right=0.75)
    par1 = host.twinx()
    host.set_xlabel('Model Complexity (%s)' % conf['complexity_label'])
    y1_label = conf['prediction_performance_label']
    y2_label = "Time (s)"
    host.set_ylabel(y1_label)
    par1.set_ylabel(y2_label)
    p1, = host.plot(complexities, mse_values, 'b-', label="prediction error")
    p2, = par1.plot(complexities, prediction_times, 'r-',
                    label="latency")
    host.legend(loc='upper right')
    host.axis["left"].label.set_color(p1.get_color())
    par1.axis["right"].label.set_color(p2.get_color())
    plt.title("Influence of varying '%s' on %s" % (conf['changing_param'],
                                                   conf['estimator'].__name__))


for conf in configurations:
    prediction_performances, prediction_times, complexities = \
        benchmark_influence(conf)
    plot_influence(conf, prediction_performances, prediction_times,
                   complexities)
plt.show()
