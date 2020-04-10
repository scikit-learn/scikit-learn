"""
===============================================================================
Comparing different split criteria for random forest regression on toy datasets
===============================================================================
An example to compare the different split criteria available for
:class:`sklearn.ensemble.RandomForestRegressor`.
Metrics used to evaluate these splitters include Mean Squared Error (MSE), a
measure of distance between the true target (`y_true`) and the predicted output
(`y_pred`), and runtime.
For visual examples of these datasets, see
:ref:`sphx_glr_auto_examples_datasets_plot_nonlinear_regression_datasets.py`.
"""

# Authors: Vivek Gopalakrishnan <vgopala4@jhu.edu>
#          Morgan Sanchez       <msanch35@jhu.edu>
# License: BSD 3 clause

import time
from itertools import product
from joblib import Parallel, delayed

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from sklearn.datasets import (
    make_independent_noise,
    make_log_regression,
    make_multiplicative_noise,
    make_sin_regression,
    make_square_regression,
)
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error

print(__doc__)

random_state = 0

###############################################################################
noise = 100.0
simulations = {
    "Logarithmic": [make_log_regression, noise],
    r"Sine Period $4\pi$": [make_sin_regression, noise],
    "Square": [make_square_regression, noise],
    "Multiplicative": [make_multiplicative_noise, None],
    "Independence": [make_independent_noise, None],
}


###############################################################################
def _train_forest(X, y, criterion):
    """Fit RandomForestRegressor with default parameters & given criterion."""
    regr = RandomForestRegressor(
        n_estimators=500, criterion=criterion, max_features="sqrt", max_depth=5
    )
    regr.fit(X, y)
    return regr


def _test_forest(X, y, regr):
    """Calculate the accuracy of the model on a heldout set."""
    y_pred = regr.predict(X)
    return mean_squared_error(y, y_pred)


def _prep_data(sim_dict, sim_name, max_n_samples, n_dimensions, n_trials):
    """Generate train and test data for all trials."""
    # Get simulation parameters and validation dataset
    sim, noise, (X_test, y_test) = simulations[sim_name]
    n_samples = int(max_n_samples)
    n_dimensions = int(n_dimensions)

    np.random.seed(random_state)
    seeds = np.random.randint(1e8, size=n_trials)

    sim_dict[sim_name] = {}
    for i in range(n_trials):
        # Sample training data
        if noise is not None:
            X_train, y_train = sim(
                n_samples=n_samples,
                n_dimensions=n_dimensions,
                noise=noise,
                random_state=seeds[i],
            )
        else:
            X_train, y_train = sim(
                n_samples=n_samples, n_dimensions=n_dimensions,
                random_state=seeds[i]
            )
        sim_dict[sim_name][i] = (
            np.copy(X_train),
            np.copy(y_train),
            np.copy(X_test),
            np.copy(y_test),
        )
    return sim_dict


###############################################################################
def main(sim_name, sim_data, n_samples, criterion, n_dimensions, n_iter):
    """Measure the performance of RandomForest under simulation conditions.
    Parameters
    ----------
    simulation_name : str
        Key from `simulations` dictionary.
    sim_data: tuple (X_train, y_train, X_test, y_test)
            X_train : array, shape (n_train_samples, n_features)
                All X training data for given simulation
            y_train : array, shape (n_train_samples, n_outputs)
                All y training data for given simulation
            X_test : array, shape (n_test_samples, n_features)
                All X testing data for given simulation
            y_test : array, shape (n_test_samples, n_outputs)
                All y testing data for given simulation
    n_samples : int
        Number of training samples.
    criterion : {'mse', 'mae', 'friedman_mse'}
        Split criterion used to train forest:
        - 'mse'
            Mean Squared Error
        - 'mae'
            Mean Absolute Error
        - 'friedman_mse'
            Friedman Mean Squared Error
    n_dimensions : int
        Number of features and targets to sample.
    n_iter : int
        Which repeat of the same simulation parameter we're on. Ignored.
    Returns
    -------
    sim_name : str
        Key from `simulations` dictionary.
    n_samples : int
        Number of training samples.
    criterion : string
        Split criterion used to train forest. Choose from
        ("mse", "mae", "friedman_mse", "axis", "oblique").
    n_dimensions : int, optional (default=10)
        Number of features and targets to sample.
    score : float
        Euclidean distance between y_pred and y_test.
    runtime : float
        Runtime (in seconds).
    """
    print(sim_name, n_samples, criterion, n_dimensions, n_iter)

    # Unpack training and testing data
    X_train, y_train, X_test, y_test = sim_data

    # Get subset of training data
    curr_X_train = X_train[0:n_samples]
    curr_y_train = y_train[0:n_samples]

    # Train forest
    start = time.time()
    regr = _train_forest(curr_X_train, curr_y_train, criterion)
    stop = time.time()

    # Evaluate on testing data and record runtime
    mse = _test_forest(X_test, y_test, regr)
    runtime = stop - start

    return (sim_name, n_samples, criterion, n_dimensions, mse, runtime)


###############################################################################
print("Constructing parameter space...")

# Declare simulation parameters
n_dimensions = 10
simulation_names = simulations.keys()
sample_sizes = np.arange(5, 51, 3)
criteria = ["mae", "mse", "friedman_mse"]

# Number of times to repeat each simulation setting
n_repeats = 30

# Create the parameter space
params = product(simulation_names, sample_sizes, criteria, range(n_repeats))


###############################################################################
print("Constructing validation datasets...")
for simulation_name, (sim, noise) in simulations.items():
    if noise is not None:
        X_test, y_test = sim(
            n_samples=1000,
            n_dimensions=n_dimensions,
            noise=noise,
            random_state=random_state,
        )
    else:
        X_test, y_test = sim(
            n_samples=1000, n_dimensions=n_dimensions,
            random_state=random_state
        )
    simulations[simulation_name].append((X_test, y_test))


###############################################################################
print("Running simulations...")

# Generate training and test data for simulations
sim_data = {}
for sim in simulation_names:
    sim_data = _prep_data(sim_data, sim, sample_sizes[-1],
                          n_dimensions, n_repeats)

# Run the simulations in parallel
data = Parallel(n_jobs=-2)(
    delayed(main)(sim_name, sim_data[sim_name][n_iter], n, crit,
                  n_dimensions, n_iter)
    for sim_name, n, crit, n_iter in params
)

# Save results as a DataFrame
columns = ["simulation", "n_samples", "criterion", "n_dimensions",
           "mse", "runtime"]
df = pd.DataFrame(data, columns=columns)

# Plot the results
sns.relplot(
    x="n_samples",
    y="mse",
    hue="criterion",
    col="simulation",
    kind="line",
    data=df,
    facet_kws={"sharey": False, "sharex": True},
)
plt.tight_layout()
plt.show()
