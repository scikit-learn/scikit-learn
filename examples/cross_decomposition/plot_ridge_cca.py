"""
==========================================
Effect of Alpha on Ridge CCA
==========================================
This example illustrates tuning Ridge CCA using a range of alpha values and
visualizes how train and test correlation and covariance vary with alpha.
"""

import matplotlib.pyplot as plt
import numpy as np

from sklearn.cross_decomposition import CCA, PLSCanonical, RidgeCCA

n = 100
p = q = 50

# generate random data with correlation
rng = np.random.default_rng(0)
z = rng.normal(size=(n, 1))
Wx = rng.normal(size=(1, p))
Wy = rng.normal(size=(1, q))

# random covariance matrix with correlated variables
Cx = rng.uniform(low=-1, high=1, size=(p, p))
Cx = np.dot(Cx, Cx.T)
Cy = rng.uniform(low=-1, high=1, size=(q, q))
Cy = np.dot(Cy, Cy.T)

X = np.dot(z, Wx) + rng.multivariate_normal(mean=np.zeros(p), cov=Cx, size=n)
Y = np.dot(z, Wy) + rng.multivariate_normal(mean=np.zeros(p), cov=Cy, size=n)

X_train = X[: n // 2]
Y_train = Y[: n // 2]
X_test = X[n // 2 :]
Y_test = Y[n // 2 :]


# Function to calculate covariance and correlation
def calculate_cov_corr(X, Y):
    cov = np.cov(X.T, Y.T)[0, 1]
    corr = np.corrcoef(X.T, Y.T)[0, 1]
    return np.mean(cov), np.mean(corr)


# Fit and evaluate CCA and PLS models
cca_model = CCA(n_components=1)
pls_model = PLSCanonical(n_components=1)

cca_model.fit(X_train, Y_train)
pls_model.fit(X_train, Y_train)

# Transform test data and calculate metrics for CCA and PLS
X_train_cca, Y_train_cca = cca_model.transform(X_train, Y_train)
X_test_cca, Y_test_cca = cca_model.transform(X_test, Y_test)
X_train_pls, Y_train_pls = pls_model.transform(X_train, Y_train)
X_test_pls, Y_test_pls = pls_model.transform(X_test, Y_test)

cca_train_cov, cca_train_corr = calculate_cov_corr(X_train_cca, Y_train_cca)
cca_test_cov, cca_test_corr = calculate_cov_corr(X_test_cca, Y_test_cca)
pls_train_cov, pls_train_corr = calculate_cov_corr(X_train_pls, Y_train_pls)
pls_test_cov, pls_test_corr = calculate_cov_corr(X_test_pls, Y_test_pls)

# Analyzing effect of alpha on train and test correlation and covariance
alphas = [1e-3, 1e-2, 1e-1, 1, 10, 100, 1000]
train_covs = []
test_covs = []
train_corrs = []
test_corrs = []

for alpha in alphas:
    model = RidgeCCA(n_components=1, alpha_x=alpha, alpha_y=alpha)
    model.fit(X_train, Y_train)

    # Train data
    X_train_trans, Y_train_trans = model.transform(X_train, Y_train)
    train_cov, train_corr = calculate_cov_corr(X_train_trans, Y_train_trans)
    train_covs.append(train_cov)
    train_corrs.append(train_corr)

    # Test data
    X_test_trans, Y_test_trans = model.transform(X_test, Y_test)
    test_cov, test_corr = calculate_cov_corr(X_test_trans, Y_test_trans)
    test_covs.append(test_cov)
    test_corrs.append(test_corr)

plt.figure(figsize=(10, 6))

# Define markers for better visibility
markers = ["o", "s", "D", "^", "v", "*", "x", "."]

# Plot with markers and different line styles
plt.plot(
    alphas,
    train_corrs,
    label="Ridge Train Correlation",
    color="blue",
    linestyle="-",
    marker=markers[0],
)
plt.plot(
    alphas,
    test_corrs,
    label="Ridge Test Correlation",
    color="blue",
    linestyle="--",
    marker=markers[1],
)
plt.axhline(cca_train_corr, color="red", linestyle="-", label="CCA Train Correlation")
plt.axhline(cca_test_corr, color="red", linestyle="--", label="CCA Test Correlation")
plt.axhline(pls_train_corr, color="green", linestyle="-", label="PLS Train Correlation")
plt.axhline(pls_test_corr, color="green", linestyle="--", label="PLS Test Correlation")

# Add gridlines
plt.grid(which="both", linestyle="--", linewidth=0.5)

# Increase font sizes
plt.xlabel("Alpha", fontsize=14)
plt.ylabel("Mean Correlation", fontsize=14)
plt.title("Effect of Alpha on Correlation", fontsize=16)

# Modify x-axis to log scale
plt.xscale("log")

# Increase tick sizes
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

# Move the legend outside of the plot
plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
# Adjust the rect to make room for the legend
plt.show()
