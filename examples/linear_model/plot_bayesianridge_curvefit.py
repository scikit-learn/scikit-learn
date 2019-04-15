import matplotlib.pyplot as plt
from sklearn.linear_model import BayesianRidge
import numpy as np
"""
============================================
Curve Fitting with Bayesian Ridge Regression
============================================

Computes a Bayesian Ridge Regression of Sinusoids.

See :ref:`bayesian_ridge_regression` for more information on the regressor.

In general, when fitting a curve with a polynomial by Bayesian ridge
regression (evidence approxymation), the selection of initial values of
hyperparameters (alpha, lambda) may be important.

In this example, the sinusoid is approximated by a polynomial using different
pairs of initial values.

When starting from the default values (alpha_init =1.90, lambda_init = 1.),
the bias of the resulting curve is large, and the variance is small.
So, lambda_init should be relatively small (1.e-3) so as to reduce the bias.

Also, by evaluating the evidences (marginal log-likelihoods, logL) of
these models, we can determine which one is better.
It can be concluded that the model with larger evidence are more likely.
"""
print(__doc__)


# #############################################################################
# Generate sinusoidal data with noise

def func(x):
    return np.sin(2*np.pi*x)


size = 25
np.random.seed(1234)
xtrain = np.random.uniform(0., 1., size)
ytrain = func(xtrain)+np.random.normal(scale=0.1, size=size)
xtest = np.linspace(0., 1., 100)


# #############################################################################
# Fit by cubic polynomial
nOrder = 3
Xtrain = np.vander(xtrain, nOrder+1, increasing=True)
Xtest = np.vander(xtest, nOrder+1, increasing=True)

# #############################################################################
# Bayesian ridge regression with different initial value pairs
inits = (1., 1.e-3)
regs = [BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True),
        BayesianRidge(tol=1e-6, fit_intercept=False, compute_score=True,
                      alpha_init=inits[0], lambda_init=inits[1])]

# #############################################################################
# Plot the true and predicted curves with logLs (evidences)
fig, ax = plt.subplots(1, 2, figsize=(8, 4))
for i, reg in enumerate(regs):
    reg.fit(Xtrain, ytrain)
    ymean, ystd = reg.predict(Xtest, return_std=True)

    ax[i].plot(xtest, func(xtest), color="blue", label="sin(2$πx$)")
    ax[i].scatter(xtrain, ytrain, s=50, alpha=0.5, label="observation")
    ax[i].plot(xtest, ymean, color="red", label="predict mean")
    ax[i].fill_between(xtest, ymean-ystd, ymean+ystd,
                       color="pink", alpha=0.5, label="predict std")
    ax[i].set_ylim(-1.3, 1.3)
    ax[i].legend()
    if i == 0:
        ax[i].set_title("$α$_init$={:.2f} ,λ$_init$={}$ (Default)".format(
            1./np.var(ytrain), 1.))
    elif i == 1:
        ax[i].set_title(
            "$α$_init$={} ,λ$_init$={}$".format(inits[0], inits[1]))
    ax[i].text(0.05, -1.0, "$α={:.1f}$\n$λ={:.3f}$\nlog$L={:.1f}$".format(
        reg.alpha_, reg.lambda_, reg.scores_[-1]), fontsize=12)

plt.tight_layout()
plt.show()
