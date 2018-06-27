"""
=================================================
Imputing missing values using multiple imputation
=================================================

By default, the ChainedImputer performs single imputation: a method where every
missing value is replaced with one imputed value. The strength of the method is
that it allows for finding unbiased statistical estimates due to its chained
character. However, the disadvantage is that every imputed value is treated as
if the value was observed, leading to an imputed dataset that does not reflect
the uncertainty that occurs due to the presence of missing values. This makes it
hard to find valid statistical inferences because the variance (and standard error)
of statistical estimates become too small.

An alternative is using the ChainedImputer to perform multiple imputation: a method
where every missing value is imputed multiple times. The procedure results in
multiple datasets where the observed data is similar in every dataset, but the imputed
data is different. All desired steps after imputation are performed on every dataset,
including the analysis. Then, Rubin's pooling rules are used to combine the estimates
into one final result.

In this example we will show how to use the ChainedImputer to perform multiple imputation,
what the effect is on the standard error of beta coefficients and how to set up a prediction
model using multiple imputation.
"""

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

from sklearn.datasets import load_boston
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.impute import SimpleImputer, ChainedImputer
from sklearn.metrics import mean_squared_error as mse

rng = np.random.RandomState(0)

def ampute(X, missing_rate = 0.75, mech = "MCAR"):

    n_samples = X.shape[0]
    n_features = X.shape[1]
    X_incomplete = X.copy()

    # MCAR mechanism
    if mech == 'MCAR':
        for i in np.arange(n_features):
            dropped_indices = np.array(np.random.choice(np.arange(n_samples), size=int(missing_rate * n_samples), replace=False))
            X_incomplete[dropped_indices[:, None], i] = None

    # MNAR mechanism
    if mech == "MNAR":
        for i in np.arange(n_features):
            data_values = -np.mean(X[:, i]) + X[:, i]
            weights = list(map(lambda x: math.exp(x) / (1 + math.exp(x)), data_values))
            probs = np.array(weights) / np.sum(np.array(weights))
            dropped_indices = np.array(np.random.choice(np.arange(n_samples), size=int(missing_rate * n_samples), p=probs, replace=False))
            X_incomplete[dropped_indices[:, None], i] = None

    return X_incomplete

def calculate_variance_of_beta_estimates(y_true, y_pred, X):

    residuals = np.sum((y_true - y_pred)**2)
    sigma_hat_squared = (1 / (len(y_true) - 2)) * residuals
    X_prime_X = np.dot(X.T, X)
    covariance_matrix = sigma_hat_squared / X_prime_X
    vars = np.diag(covariance_matrix)

    return vars

### EXAMPLE 1.
### COMPARE STATISTICAL ESTIMATES AND THEIR VARIANCE FOR LINEAR REGRESSION MODEL

def get_results_full_dataset(X, y):

    # Perform linear regression on full data as a way of comparison
    estimator = LinearRegression()
    estimator.fit(X, y)
    y_predict = estimator.predict(X)

    # Save the beta estimates
    # The variance of these estimates
    # And 1.96 * standard error of the estimates (useful to know the 95% confidence interval)
    full_coefs = estimator.coef_
    full_vars = calculate_variance_of_beta_estimates(y, y_predict, X)
    full_errorbar = 1.96 * np.sqrt(full_vars)

    return full_coefs, full_vars, full_errorbar

def get_results_chained_imputation(X_incomplete, y):

    # Impute incomplete data with ChainedImputer
    # Setting burnin at 99 and using only the last imputation
    imputer = ChainedImputer(n_burn_in=99, n_imputations=1)
    imputer.fit(X_incomplete)
    X_imputed = imputer.transform(X_incomplete)

    # Perform linear regression on chained single imputed data
    # Estimate beta estimates and their variances
    estimator = LinearRegression()
    estimator.fit(X_imputed, y)
    y_predict = estimator.predict(X_imputed)

    # Save the beta estimates
    # The variance of these estimates
    # And 1.96 * standard error of the estimates
    chained_coefs = estimator.coef_
    chained_vars = calculate_variance_of_beta_estimates(y, y_predict, X_imputed)
    chained_errorbar = 1.96 * np.sqrt(chained_vars)

    return chained_coefs, chained_vars, chained_errorbar

def get_results_mice_imputation(X_incomplete, y):

    # Impute incomplete data using the ChainedImputer as a MICEImputer
    # Setting burnin at 99, using only last imputation and loop this procedure m times
    m = 5
    multiple_imputations = []

    for i in range(m):

        imputer = ChainedImputer(n_burn_in=99, n_imputations=1,random_state=i)
        imputer.fit(X_incomplete)
        X_imputed = imputer.transform(X_incomplete)
        multiple_imputations.append(X_imputed)

    # Perform a model on each of the m imputed datasets
    # Estimate the estimates for each model/dataset
    m_coefs = []
    m_vars = []
    for i in range(m):

        estimator = LinearRegression()
        estimator.fit(multiple_imputations[i], y)
        y_predict = estimator.predict(multiple_imputations[i])

        m_coefs.append(estimator.coef_)
        m_vars.append(calculate_variance_of_beta_estimates(y, y_predict, multiple_imputations[i]))

    # Calculate the end estimates by applying Rubin's rules
    # Rubin's rules can be slightly different for different types of estimates
    # In case of linear regression, these are the rules:
    # The value of every estimate is the mean of estimates in each of the m datasets
    # The variance of these estimates is a combination of the variance of each of the m estimates (Ubar)
    # And the variance between the m estimates (B)

    Qbar = np.mean(m_coefs, axis = 0)
    Ubar = np.mean(m_vars, axis = 0)
    B = (1 / (m-1)) * np.mean((Qbar - m_coefs) ** 2, axis = 0)
    T = Ubar + B + (B/m)

    # The final 1.96 * standard error is then the sqrt of the variance
    mice_errorbar = 1.96 * np.sqrt(T)

    return Qbar, T, mice_errorbar

# The original MICE procedure includes all variables inluding the output variable in the imputation
# process. The idea is that the imputation model should at least contain the analysis model to
# result in unbiased estimates
def get_results_mice_imputation_includingy(X_incomplete, y):

    # Impute incomplete data using the ChainedImputer as a MICEImputer
    # Now using the output variable in the imputation loop
    m = 5
    multiple_imputations = []

    for i in range(m):

        Xy = np.column_stack((X_incomplete, y))
        imputer = ChainedImputer(n_burn_in=99, n_imputations=1, random_state=i)
        imputer.fit(Xy)
        data_imputed = imputer.transform(Xy)

        # We save only the X imputed data because we don't want to use y to predict y later on
        X_imputed = data_imputed[:, :-1]
        multiple_imputations.append(X_imputed)

    # Perform linear regression on mice multiple imputed data
    # Estimate beta estimates and their variances
    m_coefs = []
    m_vars = []
    for i in range(m):

        estimator = LinearRegression()
        estimator.fit(multiple_imputations[i], y)
        y_predict = estimator.predict(multiple_imputations[i])

        m_coefs.append(estimator.coef_)
        m_vars.append(calculate_variance_of_beta_estimates(y, y_predict, multiple_imputations[i]))

    # Calculate the end results by applying Rubin's rules
    # The value of every estimate is the mean of the values over the m datasets
    # The variance of these estimates is a combination of the variance of each of the m estimates (Ubar)
    # And the variance between the m estimates (B)

    Qbar = np.mean(m_coefs, axis = 0)
    Ubar = np.mean(m_vars, axis = 0)
    B = (1 / (m-1)) * np.mean((Qbar - m_coefs) ** 2, axis = 0)
    T = Ubar + B + (B/m)

    # The final 1.96 * standard error is then the sqrt of the variance
    mice_errorbar = 1.96 * np.sqrt(T)

    return Qbar, T, mice_errorbar

# Now lets run these imputation procedures
# We use the Boston dataset and analyze the outcomes of the beta coefficients and their standard errors
# We standardize the data before running the procedure to be able to compare the coefficients
# We run the procedure for 3 missingness mechanisms (MCAR, MAR and MNAR)

dataset = load_boston()
X_full, y = dataset.data, dataset.target

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X_full)
y_scaled = stats.zscore(y)

print("Executing Example 1 MCAR Missingness")
Boston_X_incomplete_MCAR = ampute(X_scaled, mech = "MCAR")

full_coefs, full_vars, full_errorbar = get_results_full_dataset(X_scaled, y_scaled)
chained_coefs, chained_vars, chained_errorbar = get_results_chained_imputation(Boston_X_incomplete_MCAR, y_scaled)
mice_coefs, mice_vars, mice_errorbar = get_results_mice_imputation(Boston_X_incomplete_MCAR, y_scaled)
mice_y_coefs, mice_y_vars, mice_y_errorbar = get_results_mice_imputation_includingy(Boston_X_incomplete_MCAR, y_scaled)

coefs = (full_coefs, chained_coefs, mice_coefs, mice_y_coefs)
vars = (full_vars, chained_vars, mice_vars, mice_y_vars)
errorbars = (full_errorbar, chained_errorbar, mice_errorbar, mice_y_errorbar)

# We plot the results
n_situations = 4
n = np.arange(n_situations)
n_labels = ['Full Data', 'Chained Imputer', 'Mice Imputer', 'Mice Imputer with y']
colors = ['r', 'orange', 'b', 'purple']
width = 0.3
plt.figure(figsize=(24, 16))

plt1 = plt.subplot(211)
for j in n:
    plt1.bar(np.arange(len(coefs[j])) + (3*j*(width/n_situations)), coefs[j], width = width, color = colors[j])
plt.legend(n_labels)

plt2 = plt.subplot(212)
for j in n:
    plt2.bar(np.arange(len(errorbars[j])) + (3*j*(width/n_situations)), errorbars[j], width = width, color = colors[j])

plt1.set_title("MCAR Missingness")
plt1.set_ylabel("Beta Coefficients")
plt2.set_ylabel("Standard Errors")
plt1.set_xlabel("Features")
plt2.set_xlabel("Features")

plt.show()

### EXAMPLE 2. ###
### SHOW MULTIPLE IMPUTATION IN PREDICTION CONTEXT ###

# In this example, we show how to apply the imputer in a train/test situation
# There are two approaches to get the end result of the prediction model
# In approach 1 you calculate the evaluation metric for every i in m and later average these values
# In approach 2 you average the predictions of every i in m and then calculate the evaluation metric

def get_results_full_data(X_train, X_test, y_train, y_test):

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)
    X_test_scaled = scaler.transform(X_test)

    estimator = LinearRegression()
    estimator.fit(X_train_scaled, y_train)
    y_predict = estimator.predict(X_test_scaled)
    mse_full = mse(y_test, y_predict)

    return mse_full

# Perform pipeline for i in m
# Approach 1: pool the mse values of the m datasets
def get_results_multiple_imputation_approach1(X_train, X_test, y_train, y_test):

    m = 5
    multiple_mses = []

    for i in range(m):

        # Fit the imputer for every i in im
        # Be aware that you fit the imputer on the train data
        # And apply to the test data
        imputer = ChainedImputer(n_burn_in=99, n_imputations=1, random_state=i)
        X_train_imputed = imputer.fit_transform(X_train)
        X_test_imputed = imputer.transform(X_test)

        # Perform the steps you wish to take before fitting the estimator
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train_imputed)
        X_test_scaled = scaler.transform(X_test_imputed)

        # Finally fit the estimator and calculate the error metric for every i in m
        estimator = LinearRegression()
        estimator.fit(X_train_scaled, y_train)
        y_predict = estimator.predict(X_test_scaled)
        mse_approach1 = mse(y_test, y_predict)
        multiple_mses.append(mse_approach1)

    # Average the error metric over the m loops to get a final result
    mse_approach1 = np.mean(multiple_mses, axis=0)

    return mse_approach1

# Approach 2: average the predictions of the m datasets and then calculate the mse
def get_results_multiple_imputation_approach2(X_train, X_test, y_train, y_test):

    m = 5
    multiple_predictions = []

    for i in range(m):

        # Fit the imputer for every i in im
        # Be aware that you fit the imputer on the train data
        # And apply to the test data
        imputer = ChainedImputer(n_burn_in=99, n_imputations=1, random_state=i)
        X_train_imputed = imputer.fit_transform(X_train)
        X_test_imputed = imputer.transform(X_test)

        # Perform the steps you wish to take before fitting the estimator
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train_imputed)
        X_test_scaled = scaler.transform(X_test_imputed)

        # Finally fit the estimator and calculate the predictions for every i in m
        estimator = LinearRegression()
        estimator.fit(X_train_scaled, y_train)
        y_predict = estimator.predict(X_test_scaled)
        multiple_predictions.append(y_predict)

    # Average the predictions over the m loops
    # Then calculate the error metric
    predictions_average = np.mean(multiple_predictions, axis=0)
    mse_approach2 = mse(y_test, predictions_average)

    return mse_approach2

def perform_simulation(dataset, X_incomplete, nsim = 10):

    X_full, y = dataset.data, dataset.target
    outcome = []

    for j in np.arange(nsim):

        train_indices, test_indices = train_test_split(np.arange(X_full.shape[0]))

        X_incomplete_train = X_incomplete[train_indices]
        X_full_train = X_full[train_indices]
        X_incomplete_test = X_incomplete[test_indices]
        X_full_test = X_full[test_indices]
        y_train = y[train_indices]
        y_test = y[test_indices]

        mse_full = get_results_full_data(X_full_train, X_full_test, y_train, y_test)
        mse_approach1 = get_results_multiple_imputation_approach1(X_incomplete_train, X_incomplete_test, y_train, y_test)
        mse_approach2 = get_results_multiple_imputation_approach2(X_incomplete_train, X_incomplete_test, y_train, y_test)

        outcome.append((mse_full, mse_approach1, mse_approach2))

    return np.mean(outcome, axis = 0), np.std(outcome, axis = 0)

# Execute
print("Executing Example 1 MCAR Missingness")
Boston_X_incomplete_MCAR = ampute(X_scaled, mech = "MCAR")
mse_means, mse_std = perform_simulation(load_boston(), Boston_X_incomplete_MCAR, nsim=10)

# Plot results
n_situations = 3
n = np.arange(n_situations)
n_labels = ['Full Data', 'Average MSE', 'Average Predictions']
colors = ['r', 'green', 'yellow']
width = 0.3
plt.figure(figsize=(6, 6))

plt1 = plt.subplot(111)
for j in n:
    plt1.bar(j, mse_means[j], yerr = mse_std[j],
             width = width, color = colors[j])

plt1.set_title("MCAR Missingness")
plt1.set_ylabel("Mean Squared Error")
plt.legend(n_labels)
plt.show()
