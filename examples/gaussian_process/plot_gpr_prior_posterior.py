"""
==========================================================================
Illustration of prior and posterior Gaussian process for different kernels
==========================================================================

This example illustrates the prior and posterior of a GPR with different
kernels. Mean, standard deviation, and 10 samples are shown for both prior
and posterior.
"""
print(__doc__)

# Authors: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
#          Guillaume Lemaitre <g.lemaitre58@gmail.com>
# License: BSD 3 clause

# %%
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (
    RBF,
    Matern,
    RationalQuadratic,
    ExpSineSquared,
    DotProduct,
    ConstantKernel,
)


# %%
import numpy as np
import matplotlib.pyplot as plt


def plot_prior_gpr(gpr_model, n_sampled_function, ax):
    """Plot the prior of the Gaussian process.

    Parameters
    ----------
    gpr_model : `GaussianProcessRegressor`
        An unfitted `GaussianProcessRegressor`.
    n_sampled_functions : int
        The number of prior functions to sample.
    ax : matplotlib axis
        The axis where to plot the prior.
    """
    # generate the data
    x = np.linspace(0, 5, 100)
    X = x.reshape(-1, 1)

    # get the priors
    y_mean, y_std = gpr_model.predict(X, return_std=True)
    y_samples = gpr_model.sample_y(X, n_sampled_function)

    for idx, single_prior in enumerate(y_samples.T):
        ax.plot(x, single_prior)
    ax.plot(x, y_mean, color="black")
    ax.fill_between(x, y_mean - y_std, y_mean + y_std, alpha=0.2, color="black")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(f"Prior (kernel:  {kernel})", fontsize=12)


def plot_posterior_gpr(gpr_model, X_train, y_train, n_sampled_function, ax):
    gpr_model.fit(X_train, y_train)

    # generate the data
    x = np.linspace(0, 5, 100)
    X = x.reshape(-1, 1)
    n_posterior = 5

    # get the posteriors
    y_mean, y_std = gp.predict(X, return_std=True)
    y_samples = gp.sample_y(X, n_posterior)

    for idx, single_prior in enumerate(y_samples.T):
        ax.plot(x, single_prior, label=f"Sampled function #{idx + 1}")
    ax.plot(x, y_mean, color="black", label="Mean prior")
    ax.fill_between(
        x,
        y_mean - y_std,
        y_mean + y_std,
        alpha=0.2,
        color="black",
        label="Mean +/- Std. Dev.",
    )
    ax.scatter(X_train[:, 0], y_train, color="red")

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(
        f"Posterior (kernel: {gpr_model.kernel})\n "
        f"Log-Likelihood: "
        f"{gpr_model.log_marginal_likelihood(gpr_model.kernel_.theta):.3f}",
        fontsize=12,
    )


# %%
rng = np.random.RandomState(4)
X_train = rng.uniform(0, 5, 10).reshape(-1, 1)
y_train = np.sin((X_train[:, 0] - 2.5) ** 2)


# %%
kernels = [
    1.0 * RBF(length_scale=1.0, length_scale_bounds=(1e-1, 10.0)),
    1.0 * RationalQuadratic(length_scale=1.0, alpha=0.1, alpha_bounds=(1e-5, 1e10)),
    1.0
    * ExpSineSquared(
        length_scale=1.0,
        periodicity=3.0,
        length_scale_bounds=(0.1, 10.0),
        periodicity_bounds=(1.0, 10.0),
    ),
    ConstantKernel(0.1, (0.01, 10.0))
    * (DotProduct(sigma_0=1.0, sigma_0_bounds=(0.1, 10.0)) ** 2),
    1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0), nu=1.5),
]

# %%
_, axs = plt.subplots(nrows=2, figsize=(6, 8))
gp = GaussianProcessRegressor(kernel=kernels[0])
plot_prior_gpr(gp, ax=axs[0])
plot_posterior_gpr(gp, X_train, y_train, ax=axs[1])
axs[1].legend(bbox_to_anchor=(1.05, 1.5), loc="upper left")
plt.subplots_adjust(hspace=0.4)

# %%
kernels = [
    1.0 * RBF(length_scale=1.0, length_scale_bounds=(1e-1, 10.0)),
    1.0 * RationalQuadratic(length_scale=1.0, alpha=0.1, alpha_bounds=(1e-5, 1e10)),
    1.0
    * ExpSineSquared(
        length_scale=1.0,
        periodicity=3.0,
        length_scale_bounds=(0.1, 10.0),
        periodicity_bounds=(1.0, 10.0),
    ),
    ConstantKernel(0.1, (0.01, 10.0))
    * (DotProduct(sigma_0=1.0, sigma_0_bounds=(0.1, 10.0)) ** 2),
    1.0 * Matern(length_scale=1.0, length_scale_bounds=(1e-1, 10.0), nu=1.5),
]

for kernel in kernels:
    # Specify Gaussian Process
    gp = GaussianProcessRegressor(kernel=kernel)

    # Plot prior
    plt.figure(figsize=(8, 8))
    plt.subplot(2, 1, 1)
    X_ = np.linspace(0, 5, 100)
    y_mean, y_std = gp.predict(X_[:, np.newaxis], return_std=True)
    plt.plot(X_, y_mean, "k", lw=3, zorder=9)
    plt.fill_between(X_, y_mean - y_std, y_mean + y_std, alpha=0.2, color="k")
    y_samples = gp.sample_y(X_[:, np.newaxis], 10)
    plt.plot(X_, y_samples, lw=1)
    plt.xlim(0, 5)
    plt.ylim(-3, 3)
    plt.title("Prior (kernel:  %s)" % kernel, fontsize=12)

    # Generate data and fit GP
    rng = np.random.RandomState(4)
    X = rng.uniform(0, 5, 10)[:, np.newaxis]
    y = np.sin((X[:, 0] - 2.5) ** 2)
    gp.fit(X, y)

    # Plot posterior
    plt.subplot(2, 1, 2)
    X_ = np.linspace(0, 5, 100)
    y_mean, y_std = gp.predict(X_[:, np.newaxis], return_std=True)
    plt.plot(X_, y_mean, "k", lw=3, zorder=9)
    plt.fill_between(X_, y_mean - y_std, y_mean + y_std, alpha=0.2, color="k")

    y_samples = gp.sample_y(X_[:, np.newaxis], 10)
    plt.plot(X_, y_samples, lw=1)
    plt.scatter(X[:, 0], y, c="r", s=50, zorder=10, edgecolors=(0, 0, 0))
    plt.xlim(0, 5)
    plt.ylim(-3, 3)
    plt.title(
        "Posterior (kernel: %s)\n Log-Likelihood: %.3f"
        % (gp.kernel_, gp.log_marginal_likelihood(gp.kernel_.theta)),
        fontsize=12,
    )
    plt.tight_layout()

plt.show()

# %%
