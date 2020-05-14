"""
===============================================================================
Comparing Different Split Criteria for Random Forest Regression on Toy Datasets
===============================================================================
An example to compare the different split criteria available for
:class:`sklearn.ensemble.RandomForestRegressor`.
Metrics used to evaluate these splitters include runtime and Mean Squared Error (MSE), a
measure of distance between the true target (`y_true`) and the predicted output
(`y_pred`).  
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

###############################################################################
# Initialize Random State, Noise, and Simulation Dictionary
# ----------------------------------------------------------
# 
# For each simulation, we will set the appropriate data generation function and 
# noise amount. Note that Multiplicative and Independence data will be noiseless.
# Please see :ref:`sphx_glr_auto_examples_datasets_plot_nonlinear_regression_datasets.py`
# for more information on the datasets used in this example.

random_state = 0
noise = 100.0
simulations = {
    "Logarithmic": [make_log_regression, noise],
    r"Sine Period $4\pi$": [make_sin_regression, noise],
    "Square": [make_square_regression, noise],
    "Multiplicative": [make_multiplicative_noise, None],
    "Independence": [make_independent_noise, None],
}

###############################################################################
# Define a Data Preparation Function
# -----------------------------------
#
# This function generates train and test data for all trials. 
#
# First, we obtain an array of seeds for each trial.
#
# Then, for the provided simulation name (i.e. "Logarithmic", "Square", etc),
# we create a dictionary containing X_train and y_train for each trial
# using that trial's random seed in addition to a single X_test and y_test
# that is consistent across all trials. To obtain this data we use
# the simulation name's associated make_* function from the above section.

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
# Define a Function to Train a Forest
# -----------------------------------
#
# Given X and y (training data) and a split criterion, this function fits a 
# random forest regressor with 500 trees, max_features of sqrt and max_depth of 5.

def _train_forest(X, y, criterion):
    """Fit RandomForestRegressor with selected parameters & given criterion."""
    regr = RandomForestRegressor(
        n_estimators=500, criterion=criterion, max_features="sqrt", max_depth=5
    )
    regr.fit(X, y)
    return regr

###############################################################################
# Define a Function to Test a Forest
# ----------------------------------
#
# Given X_test, y_test (which are the true predictions), and regr (the random forest
# regressor), this function returns the Mean Squared Error. In other words, it assesses
# the performance of the forest.

def _test_forest(X, y, regr):
    """Calculate the accuracy of the model on a heldout set."""
    y_pred = regr.predict(X)
    return mean_squared_error(y, y_pred)

###############################################################################
# Define a Function That Performs a Single Trial
# ----------------------------------------------
#
# For a given simulation name, trial, number of samples, and criterion, 
# this function trains and evaluates a random forest regressor based on MSE
# and runtime.

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

    # Unpack training and testing data
    X_train, y_train, X_test, y_test = sim_data

    # Get subset of training data
    curr_X_train = X_train[0:n_samples]
    curr_y_train = y_train[0:n_samples]

    # Train forest
    start = time.process_time()
    regr = _train_forest(curr_X_train, curr_y_train, criterion)
    stop = time.process_time()

    # Evaluate on testing data and record runtime
    mse = _test_forest(X_test, y_test, regr)
    runtime = stop - start

    return (sim_name, n_samples, criterion, n_dimensions, mse, runtime)


###############################################################################
# Constructing the Parameter Space
# ---------------------------------
#
# For this experiment, we will be training a random forest for each simulation
# type, sample size, and split criterion 8 times. 
#
# The targets will be 10 dimensional. In other words, there will be 10 predictors
# because the aim of this experiment is to provide a nonlinear multioutput regression
# example.

print("Constructing parameter space...")

# Declare simulation parameters
n_dimensions = 10
simulation_names = simulations.keys()
sample_sizes = np.arange(5, 51, 5)
criteria = ["mae", "mse", "friedman_mse"]

# Number of times to repeat each simulation setting
n_repeats = 8

# Create the parameter space
params = product(simulation_names, sample_sizes, criteria, range(n_repeats))


###############################################################################
# Construct Validation Datasets
# ------------------------------
#
# Here, we generate X_test and y_test. For each simulation name, we generate 
# 1000 samples each with 10 targets.

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
# Run the Experiment
# ------------------
# For each simulation type (Logarithmic, Sine, Square, Multiplicative, Independence)
# 
# - And For each split criterion (mse, mae, friedman mse)
# 
#   - Generate 8 noisy training sets (10 dimensions, 50 samples)
# 
#   - Generate 1 noisy test set (10 dimensions, 1000 samples)
# 
#   - Train and evaluate using mse and runtime for all 8 training sets with random 
#   forests (500 trees) varying number of samples ( `np.arange(5, 51, 3)` )

print("Running simulations...")

# Generate training and test data for simulations
sim_data = {}
for sim in simulation_names:
    sim_data = _prep_data(sim_data, sim, sample_sizes[-1],
                          n_dimensions, n_repeats)

# Run the simulations in parallel
data = Parallel(n_jobs=-1)(
    delayed(main)(sim_name, sim_data[sim_name][n_iter], n, crit,
                  n_dimensions, n_iter)
    for sim_name, n, crit, n_iter in params
)

# Save results as a DataFrame
columns = ["Simulation", "n_samples", "Criterion", "n_dimensions",
           "mse", "Runtime"]
df = pd.DataFrame(data, columns=columns)

###############################################################################
# Plot the Results
# ----------------
#
# The first plot shows the MSE comparison across split criteria. 
# We use this plot to determine how increasing the number of training
# samples affects the prediction error. 

sns.relplot(
    x="n_samples",
    y="mse",
    hue="Criterion",
    col="Simulation",
    kind="line",
    data=df,
    facet_kws={"sharey": False, "sharex": True},
)
plt.tight_layout()
plt.show()

###############################################################################
# The second plot shows how increasing the number of training samples affects the 
# runtime, or more specifically, the process time. 

sns.relplot(
    x="n_samples",
    y="Runtime",
    hue="Criterion",
    col="Simulation",
    kind="line",
    data=df,
    facet_kws={"sharey": False, "sharex": True},
)
plt.tight_layout()
plt.show()

###############################################################################
# Analyze the Plots
# ------------------
#
# As can be seen in the first plot, MAE seems to slightly outperform MSE and 
# Friedman_MSE split criteria in terms of prediction error for almost all
# simulation types. This slight difference in performance, however, is only on average
# and in some cases, the other split criteria outperform it.
#
# The second plot clearly shows that the runtime greatly increases with number of 
# samples for the MAE split criterion however. MSE and Friedman_MSE, on the other hand,
# are not as affected by sample size. In other words, if you care at all about 
# computation time, you may want to choose MSE or Friedman_MSE over MAE because 
# the prediction error decrease that MAE may (or may not) give you most likely 
# won't outweigh the computational cost.
#
# ---------
# Takeaways
# ---------
#
# For the most part, if your data resembles these nonlinear types, you may
# want to consider using Friedman_MSE or MSE for your split criterion because
# the performance benefits may not outweight the computational cost. Specifically,
# for independence or square data, the choice is more clear since the average MAE
# mean squared error is not better than MSE and Friedman_MSE. For multiplicative 
# data, however, the performance enhancement may be enough to sacrifice the resources
# for, but ultimately, that choice is up to you.
#
# This example provides a simple method for comparing split criteria and ultimately
# choosing the one that best fits your needs.