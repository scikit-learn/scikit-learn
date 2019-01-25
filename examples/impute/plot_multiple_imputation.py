"""
=================================================
Imputing missing values using multiple imputation
=================================================

The :class:`sklearn.impute.IterativeImputer` is able to generate multiple
imputations of the same incomplete dataset. We can then learn a regression
or classification model on different imputations of the same dataset.
This allows us to quantify how much the regressor varies due to the uncertainty
inherent in missing values. Specifically, we can estimate the variance of
predictions, or of the model coefficients, using multiple imputation.

Setting :class:`sklearn.impute.IterativeImputer`'s `sample_posterior=True` will
randomly draw values to fill each missing value from the Gaussian posterior of
:class:`sklearn.linear_model.BayesianRidge` predictions. If each
:class:`sklearn.impute.IterativeImputer` uses a different `random_state`, this
results in multiple imputations, each of which can be used to train a
predictive model.

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
"""
print(__doc__)

# Authors: Rianne Schouten <r.m.schouten@uu.nl>
#          Sergey Feldman <sergeyfeldman@gmail.com>
# License: BSD 3 clause

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats, special

from sklearn.datasets import load_boston
from sklearn.linear_model import BayesianRidge
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.impute import IterativeImputer
from sklearn.metrics import mean_squared_error as mse

rng = np.random.RandomState(0)

N_ITER = 20  # number of iterations in IterativeImputer
N_SIM = 3  # number of simulations in Example 2
ESTIMATOR = Pipeline(steps=[('scaler', StandardScaler()),
                            ('regressor', BayesianRidge())])
ESTIMATOR_IMPUTER = Pipeline(steps=[('imputer', IterativeImputer()),
                                    ('scaler', StandardScaler()),
                                    ('regressor', BayesianRidge())])


def ampute(X, missing_rate=0.75, mech="MCAR"):
    """Insert missing data into X.

    Ampute is the inverse of impute and serves at simulating a dataset with
    missing data. Two strategies are implemented to remove data: (1) missing
    completely at random (MCAR) and (2) missing not completely at random
    (MNCAR).

    Parameters
    ----------
    X : ndarray, shape (n_samples, n_features)
        The data without missing samples.

    missing_rate : float, default=0.75
        The amount of missing data in ``X``.

    mech : str, default='MCAR'
        The strategy for removing data. Must be 'MCAR' or 'NMCAR'.
        For NMCAR: values that are farther from the feature mean are more
        likely to be made missing.

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


def variance_model_coef(y_true, y_pred, X):
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
    y_predict = ESTIMATOR.fit(X, y).predict(X)

    # Save the beta estimates, the variance of these estimates and 1.96 *
    # standard error of the estimates. The latter is useful to know the 95%
    # confidence interval.
    full_coefs = ESTIMATOR.named_steps['regressor'].coef_
    full_vars = variance_model_coef(y, y_predict, X)

    return full_coefs, full_vars


def get_results_chained_imputation(X_incomplete, y, random_state=0,
                                   impute_with_y=False):
    # Impute incomplete data with IterativeImputer using single imputation
    # We perform N_ITER imputations and only use the last imputation.
    imputer = IterativeImputer(n_iter=N_ITER,
                               sample_posterior=True,
                               random_state=random_state)
    if impute_with_y:
        Xy = np.column_stack((X_incomplete, y_scaled))
        # impute Xy, but exclude last column for subsequent regression
        X_imputed = imputer.fit_transform(Xy)[:, :-1]
    else:
        X_imputed = imputer.fit_transform(X_incomplete)

    # Perform linear regression on chained single imputed data
    # Estimate beta estimates and their variances
    y_predict = ESTIMATOR.fit(X_imputed, y).predict(X_imputed)

    # Save the beta estimates, the variance of these estimates
    chained_coefs = ESTIMATOR.named_steps['regressor'].coef_
    chained_vars = variance_model_coef(y, y_predict, X_imputed)

    return chained_coefs, chained_vars


def coef_var_mice_imputation(X_incomplete, y, n_imputations=5,
                             impute_with_y=False):
    # Train a model on each of the `m` imputed datasets
    # Estimate the estimates for each model/dataset
    m_coefs = []
    m_vars = []
    for i in range(n_imputations):
        m_coef, m_var = get_results_chained_imputation(
            X_incomplete, y, random_state=i, impute_with_y=impute_with_y)
        m_coefs.append(m_coef)
        m_vars.append(m_var)

    # Calculate the end estimates by applying Rubin's rules.
    Qbar, T = rubins_pooling_rules(m_coefs, m_vars)

    return Qbar, T


###############################################################################
# Now let's run all these imputation procedures.
# We use the Boston dataset and analyze the outcomes of the beta coefficients
# and their standard errors. We standardize the data before running the
# procedure to be able to compare the coefficients. We run the procedure for
# MCAR missingness only.
#
# Note: the original multiple imputation procedure as developed under the name
# MICE includes all variables in the imputation process; including the output
# variable. The reason to do this is that the imputation model should at least
# contain the analysis model to result in unbiased estimates. In this function,
# we will also include `y` in the imputation process.


# Loading the data
X_full, y = load_boston(return_X_y=True)

# Standardizing the data
y_scaled = stats.zscore(y)

# Start the procedure
print("Executing Example 1 MCAR Missingness...")

# First, make the data incomplete with a MCAR mechanism.
X_incomplete = ampute(X_full, mech="MCAR")

# Second, run all the imputation procedures as described above.
full_coefs, full_vars = get_results_full_dataset(X_full, y_scaled)
chained_coefs, chained_vars = get_results_chained_imputation(X_incomplete,
                                                             y_scaled)
mice_coefs, mice_vars = coef_var_mice_imputation(X_incomplete, y_scaled)
mice_y_coefs, mice_y_vars = coef_var_mice_imputation(X_incomplete, y_scaled,
                                                     impute_with_y=True)

# Combine the results from the four imputation procedures.
coefs = [full_coefs, chained_coefs, mice_coefs, mice_y_coefs]
standard_errors = [1.96 * np.sqrt(v) for v in
                   (full_vars, chained_vars, mice_vars, mice_y_vars)]

# And plot the results
n_situations = 4
n = np.arange(n_situations)
n_labels = ['Full Data', 'IterativeImputer',
            'Mice Imputer', 'Mice Imputer with y']
colors = ['r', 'orange', 'b', 'purple']
width = 0.3
plt.figure(figsize=(6, 8))

plt1 = plt.subplot(211)
for j in n:
    plt1.bar(np.arange(len(coefs[j])) + (3*j*(width/n_situations)),
             coefs[j], width=width, color=colors[j])
plt.legend(n_labels)

plt2 = plt.subplot(212)
for j in n:
    plt2.bar(np.arange(len(standard_errors[j])) + (3*j*(width/n_situations)),
             standard_errors[j], width=width, color=colors[j])

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


# Use the IterativeImputer as a single imputation procedure.
def get_mse_single_imputation(X_train, X_test, y_train, y_test,
                              random_state=0):
    ESTIMATOR_IMPUTER.set_params(imputer__n_iter=N_ITER,
                                 imputer__sample_posterior=True,
                                 imputer__random_state=random_state)
    y_predict = ESTIMATOR_IMPUTER.fit(X_train, y_train).predict(X_test)

    return mse(y_test, y_predict), y_predict


# We average the predictions of the m datasets.
def get_mse_multiple_imputation_approach(X_train, X_test, y_train, y_test):
    num_mice_runs = 5
    multiple_predictions = []
    for i in range(num_mice_runs):
        _, y_predict = get_mse_single_imputation(
            X_train, X_test, y_train, y_test, random_state=i)
        multiple_predictions.append(y_predict)

    # Average the predictions over the m loops
    # Then calculate the error metric.
    predictions_average = np.mean(multiple_predictions, axis=0)
    mse_multiple = mse(y_test, predictions_average)

    return mse_multiple, predictions_average


def perform_simulation(X_full, y, X_incomplete, n_sim=10):
    # Start a simulation process that executes the process n_sim times.
    outcome = []
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
        # Apply the regression model on the full dataset as a way of comparison
        y_predict = ESTIMATOR.fit(X_full_train, y_train).predict(X_full_test)
        mse_full = mse(y_test, y_predict)
        mse_single, _ = get_mse_single_imputation(
                X_incomplete_train, X_incomplete_test, y_train, y_test)
        mse_multiple, _ = get_mse_multiple_imputation_approach(
                X_incomplete_train, X_incomplete_test, y_train, y_test)

        # Save the outcome of every simulation round
        outcome.append((mse_full, mse_single, mse_multiple))

    # Return the mean and standard deviation of the n_sim outcome values
    return np.mean(outcome, axis=0), np.std(outcome, axis=0)


# Execute the simulation
print("Executing Example 2 MCAR Missingness...")

# Perform the simulation
mse_means, mse_std = perform_simulation(X_full, y, X_incomplete, n_sim=N_SIM)

# Plot results
n_situations = 3
n = np.arange(n_situations)
n_labels = ['Full Data', 'Single Imputation', 'Multiple Imputations']
colors = ['r', 'orange', 'green']

plt.figure(figsize=(12, 6))
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
