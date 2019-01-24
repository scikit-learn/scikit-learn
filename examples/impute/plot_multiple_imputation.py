"""
=================================================
Imputing missing values using multiple imputation
=================================================

By default, :class:`sklearn.impute.IterativeImputer` performs single
imputation: a method where every missing value is replaced with one imputed
value. The chained character of the method and the possibility to draw
imputation values from the posterior distribution of a Bayesian imputation
model allows for the finding of unbiased statistical estimates. However, the
disadvantage is that every imputed value is treated as if the value was
observed, leading to an imputed dataset that does not reflect the uncertainty
that occurs due to the presence of missing values. This makes it hard to find
valid statistical inferences because the variance (and standard error) of
statistical estimates become too small.

An alternative is using :class:`sklearn.impute.IterativeImputer` to perform
multiple imputation: a method where every missing value is imputed multiple
times. The procedure results in multiple datasets where the observed data is
identical in every dataset, but the imputed data is different.
All desired steps after imputation are performed on every dataset, such as
standardization, feature engineering and model fitting.

One final model is obtained by combining the estimates of each model with
Rubin's pooling rules, which are as follows. The overall point estimate after
multiple imputation (denoted by Qbar) is the average of all the m point
estimates. The variance of the overall point estimate is a combination of
so-called within imputation variance (Ubar) and between imputation
variance (B). Ubar is the average of the m variances of the m point estimates.
Both Qbar and Ubar are corrected with a factor 1 / m to account for sampling
variance. The between imputation variance (B) is the sum of the squared
difference between Qbar and the m point estimates, corrected with a factor
1 / (m – 1). Then, the total variance (T) of the MI overall point estimate is
Ubar + B + B/m.

These rules assume that the parameters of interest are normally distributed
which is the case with, for example, estimates of the mean and regression
coefficients. Other parameters, such as correlation coefficients need
transformation to suit the assumption of normality. If it is not possible to
approximate a normal distribution, it is better to use robust summary measures
such as medians or ranges instead of using Rubin's pooling rules. This applies
to an estimate like explained variance.

In this Example we show how to use :class:`sklearn.impute.IterativeImputer`
to perform multiple imputation. In Example 1 we show the effect of Rubin’s
pooling rules on the variance of regression estimates. Due to the between
imputation variance, the standard errors of all regression coefficients are
larger with multiple imputation than with single imputation.

In Example 2 we show how to set up a prediction model using multiple
imputation. We combine the predictions from multiple imputations and show that
this results in better regression performance.

# Authors: Rianne Schouten <r.m.schouten@uu.nl>
           Sergey Feldman <sergeyfeldman@gmail.com>
# License: BSD 3 clause
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, special

from sklearn.datasets import load_boston
from sklearn.linear_model import RidgeCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.impute import IterativeImputer
from sklearn.metrics import mean_squared_error as mse

rng = np.random.RandomState(0)

N_ITER = 20  # number of iterations in IterativeImputer
N_SIM = 3  # number of simulations in Example 2


def ampute(X, missing_rate=0.75, mech="MCAR"):
    """Insert missing data into X.

    Ampute is the inverse of impute and serves at simulating a dataset with
    missing data. Two strategies are implemented to remove data: (1) missing
    completely at random (MCAR) and (2) missing not completely at random
    (MNCAR).

    For NMCAR: values that are farther from the feature mean are more likely to
    be made missing.

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features)
        The data without missing samples.

    missing_rate : float, default=0.75
        The amount of missing data in ``X``.

    mech : str, default='MCAR'
        The strategy for removing data. Must be 'MCAR' or 'NMCAR'.

    Returns
    -------
    X_out : ndarray, shape (n_samples, n_features)
        Output data with missing entries.
    """
    X_out = X.copy()
    n_drop_per_feature = int(missing_rate * X.shape[0])
    for x in np.transpose(X_out):
        # insert missing values for each feature
        if mech == 'MCAR':
            prob = None
        elif mech == 'NMCAR':
            weights = special.expit(np.abs(x - x.mean()))
            prob = weights / weights.sum()
        drop_idx = np.random.choice(X.shape[0], p=prob, replace=False,
                                    size=n_drop_per_feature)
        x[drop_idx] = np.nan
    return X_out


def calculate_variance_of_beta_estimates(y_true, y_pred, X):
    """Calculates variance of the beta estimates.
    """
    sum_sq_errors = np.sum((y_true - y_pred)**2)
    sigma_hat_squared = sum_sq_errors / (len(y_true) - 2)
    covariance_matrix = sigma_hat_squared / np.dot(X.T, X)
    return np.diag(covariance_matrix)


def rubins_pooling_rules(m_estimates, m_variances):
    """Applies Rubin's pooling rules.

    The value of every estimate is the mean of the estimates in each of the m
    datasets (Qbar). The variance of these estimates is a combination of the
    variance of each of the m estimates (Ubar) and the variance between the m
    estimates (B).
    """
    m = len(m_estimates)
    Qbar = np.sum(m_estimates, axis=0) / m
    Ubar = np.sum(m_variances, axis=0) / m
    B = np.sum((Qbar - m_estimates) ** 2, axis=0) / (m - 1)
    T = Ubar + B + (B/m)
    return Qbar, T


###############################################################################
# EXAMPLE 1. COMPARE STATISTICAL ESTIMATES AND THEIR VARIANCE USING MULTIPLE
# IMPUTATION IN A LINEAR REGRESSION MODEL.


def get_results_full_dataset(X, y):
    # Perform linear regression on full data as a way of comparison
    estimator = RidgeCV()
    estimator.fit(X, y)
    y_predict = estimator.predict(X)

    # Save the beta estimates, the variance of these estimates and 1.96 *
    # standard error of the estimates. The latter is useful to know the 95%
    # confidence interval.
    full_coefs = estimator.coef_
    full_vars = calculate_variance_of_beta_estimates(y, y_predict, X)

    return full_coefs, full_vars


def get_results_chained_imputation(X_incomplete, y, random_state=0):
    # Impute incomplete data with IterativeImputer using single imputation
    # We perform N_ITER imputations and only use the last imputation.
    imputer = IterativeImputer(n_iter=N_ITER,
                               sample_posterior=True,
                               random_state=random_state)
    imputer.fit(X_incomplete)
    X_imputed = imputer.transform(X_incomplete)

    # Perform linear regression on chained single imputed data
    # Estimate beta estimates and their variances
    estimator = RidgeCV()
    estimator.fit(X_imputed, y)
    y_predict = estimator.predict(X_imputed)

    # Save the beta estimates, the variance of these estimates
    chained_coefs = estimator.coef_
    chained_vars = calculate_variance_of_beta_estimates(
            y, y_predict, X_imputed)

    return chained_coefs, chained_vars


def get_results_mice_imputation(X_incomplete, y, num_imputations=5):
    # Perform a model on each of the m imputed datasets
    # Estimate the estimates for each model/dataset
    m_coefs = []
    m_vars = []
    for i in range(num_imputations):
        m_coef, m_var = get_results_chained_imputation(X_incomplete, y,
                                                       random_state=i)
        m_coefs.append(m_coef)
        m_vars.append(m_var)

    # Calculate the end estimates by applying Rubin's rules.
    Qbar, T = rubins_pooling_rules(m_coefs, m_vars)

    return Qbar, T


# The original multiple imputation procedure as developed under the name
# MICE includes all variables in the imputation process; including the output
# variable. The reason to do this is that the imputation model should at least
# contain the analysis model to result in unbiased estimates. In this function,
# we will also include y in the imputation process.
def get_results_mice_imputation_includingy(X_incomplete, y, num_imputations=5):
    Xy = np.column_stack((X_incomplete, y))
    Qbar, T = get_results_mice_imputation(Xy, y, num_imputations)
    return Qbar, T


# Now lets run all these imputation procedures.
# We use the Boston dataset and analyze the outcomes of the beta coefficients
# and their standard errors. We standardize the data before running the
# procedure to be able to compare the coefficients. We run the procedure for
# MCAR missingness only.
#
# Loading the data
dataset = load_boston()
X_full, y = dataset.data, dataset.target

# Standardizing the data
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_full)
y_scaled = stats.zscore(y)

# Start the procedure
print("Executing Example 1 MCAR Missingness...")

# First, make the data incomplete with a MCAR mechanism.
X_incomplete = ampute(X_scaled, mech="MCAR")

# Second, run all the imputation procedures as described above.
full_coefs, full_vars = get_results_full_dataset(X_scaled, y_scaled)
chained_coefs, chained_vars = get_results_chained_imputation(
        X_incomplete, y_scaled)
mice_coefs, mice_vars = get_results_mice_imputation(
        X_incomplete, y_scaled)
mice_y_coefs, mice_y_vars = get_results_mice_imputation_includingy(
        X_incomplete, y_scaled)

# Combine the results from the four imputation procedures.
coefs = [full_coefs, chained_coefs, mice_coefs, mice_y_coefs]
vars = [full_vars, chained_vars, mice_vars, mice_y_vars]
errorbars = [1.96 * np.sqrt(v) for v in vars]

# And plot the results
n_situations = 4
n = np.arange(n_situations)
n_labels = ['Full Data', 'IterativeImputer',
            'Mice Imputer', 'Mice Imputer with y']
colors = ['r', 'orange', 'b', 'purple']
width = 0.3
plt.figure(figsize=(24, 32))

plt1 = plt.subplot(211)
for j in n:
    plt1.bar(np.arange(len(coefs[j])) + (3*j*(width/n_situations)),
             coefs[j], width=width, color=colors[j])
plt.legend(n_labels)

plt2 = plt.subplot(212)
for j in n:
    plt2.bar(np.arange(len(errorbars[j])) + (3*j*(width/n_situations)),
             errorbars[j], width=width, color=colors[j])

plt1.set_title("MCAR Missingness")
plt1.set_ylabel("Beta Coefficients")
plt2.set_ylabel("Standard Errors")
plt1.set_xlabel("Features")
plt2.set_xlabel("Features")
plt.show()


###############################################################################
# EXAMPLE 2. SHOW MULTIPLE IMPUTATION IN A PREDICTION CONTEXT.


# In this example, we show how to apply multiple imputation in a train/test
# situation. There are two approaches to get the end result of the prediction
# model. In approach 1 you calculate the evaluation metric for every i in m and
# later average these values. In approach 2 you average the predictions of
# every i in m and then calculate the evaluation metric. We test both
# approaches.

# Apply the regression model on the full dataset as a way of comparison.
def get_mse(X_train, X_test, y_train, y_test):
    # Standardize data
    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    # Perform estimation and prediction
    estimator = RidgeCV()
    estimator.fit(X_train_scaled, y_train)
    y_predict = estimator.predict(X_test_scaled)

    return mse(y_test, y_predict), y_predict


# Use the IterativeImputer as a single imputation procedure.
def get_mse_single_imputation(X_train, X_test, y_train, y_test):
    imputer = IterativeImputer(n_iter=N_ITER,
                               sample_posterior=True,
                               random_state=0)
    X_train_imputed = imputer.fit_transform(X_train)
    X_test_imputed = imputer.transform(X_test)

    mse_single, _ = get_mse(X_train_imputed, X_test_imputed, y_train, y_test)

    return mse_single


# We average the predictions of the m datasets.
def get_mse_multiple_imputation_approach(X_train, X_test, y_train, y_test):
    m = 5
    multiple_predictions = []
    for i in range(m):
        # Fit the imputer for every i in m
        # Be aware that you fit the imputer on the train data
        # And apply to the test data
        imputer = IterativeImputer(n_iter=N_ITER,
                                   sample_posterior=True,
                                   random_state=i)
        X_train_imputed = imputer.fit_transform(X_train)
        X_test_imputed = imputer.transform(X_test)

        _, y_predict = get_mse(X_train_imputed, X_test_imputed, y_train,
                               y_test)
        multiple_predictions.append(y_predict)

    # Average the predictions over the m loops
    # Then calculate the error metric.
    predictions_average = np.mean(multiple_predictions, axis=0)
    mse_multiple = mse(y_test, predictions_average)

    return mse_multiple


def perform_simulation(dataset, X_incomplete, n_sim=10):
    X_full, y = dataset.data, dataset.target
    outcome = []

    # Start a simulation process that executes the process n_sim times.
    for j in range(n_sim):
        # First, split the data in train and test dataset.
        train_indices, test_indices = train_test_split(
                np.arange(X_full.shape[0]), random_state=j)
        X_incomplete_train = X_incomplete[train_indices]
        X_full_train = X_full[train_indices]
        X_incomplete_test = X_incomplete[test_indices]
        X_full_test = X_full[test_indices]
        y_train = y[train_indices]
        y_test = y[test_indices]

        # Second, perform the imputation procedures and calculation of the
        # error metric for every one of the  situations.
        mse_full = get_mse(X_full_train, X_full_test, y_train, y_test)
        mse_single = get_mse_single_imputation(
                X_incomplete_train, X_incomplete_test, y_train, y_test)
        mse_multiple = get_mse_multiple_imputation_approach(
                X_incomplete_train, X_incomplete_test, y_train, y_test)

        # Save the outcome of every simulation round
        outcome.append((mse_full[0], mse_single, mse_multiple))

    # Return the mean and standard deviation of the n_sim outcome values
    return np.mean(outcome, axis=0), np.std(outcome, axis=0)


# Execute the simulation
print("Executing Example 2 MCAR Missingness...")

# Perform the simulation
mse_means, mse_std = perform_simulation(dataset, X_incomplete, n_sim=N_SIM)

# Plot results
n_situations = 3
n = np.arange(n_situations)
n_labels = ['Full Data', 'Single Imputation', 'Multiple Imputations']
colors = ['r', 'orange', 'green']

plt.figure(figsize=(24, 12))
ax1 = plt.subplot(111)
for j in n:
    ax1.barh(j, mse_means[j], xerr=mse_std[j],
             color=colors[j], alpha=0.6, align='center')

ax1.set_title('MCAR Missingness')
ax1.set_yticks(n)
ax1.set_xlabel('Mean Squared Error')
ax1.invert_yaxis()
ax1.set_yticklabels(n_labels)
plt.show()
