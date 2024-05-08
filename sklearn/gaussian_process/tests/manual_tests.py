import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

from sklearn.gaussian_process import TProcessRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import (
    RBF,
    DotProduct,
    ExpSineSquared,
    WhiteKernel,
)
from sklearn.gaussian_process.kernels import (
    ConstantKernel as C,
)


def normConf(std):
    """ Returns the 95% confidence interval of the normal distribution """
    zScore = norm(0, 1).ppf(0.975)
    confidence = zScore * np.sqrt(std)
    return confidence


vecNormConf = np.vectorize(normConf)


def target_generator(X, add_noise=False):
    target = 0.5 + np.sin(3 * X)
    if add_noise:
        rng = np.random.RandomState(1)
        target += rng.normal(0, 0.3, size=target.shape)
    return target.squeeze()


def m_test():
    ### GP ###
    kernel = 1.0 * RBF(length_scale=1e-1, length_scale_bounds=(1e-2, 1e3)) + WhiteKernel(
        noise_level=1e-2, noise_level_bounds=(1e-10, 1e1)
    )
    gpr = GaussianProcessRegressor(kernel=kernel, alpha=0.0, optimizer=None)
    tpr = TProcessRegressor(kernel=kernel, v=5, alpha=0.0, optimizer=None)
    plotSPR2(gpr)


def plotSPR1(spr: GaussianProcessRegressor):
    ### Plotting Space ###
    X = np.linspace(0, 5, num=100).reshape(-1, 1)
    y = target_generator(X, add_noise=False)

    rng = np.random.RandomState(0)
    X_train = rng.uniform(0, 5, size=20).reshape(-1, 1)
    y_train = target_generator(X_train, add_noise=True)

    ### Predict Using Stochastic Process ###
    spr.fit(X_train, y_train)
    y_mean, y_std = spr.predict(X, return_std=True)

    ### Ground Truth and Observations ###
    plt.plot(X, y, label="Ground Truth", color='red')
    plt.scatter(
        x=X_train[:, 0],
        y=y_train,
        color="black",
        alpha=0.4,
        label="Observations",
    )

    ### Stochastic Process ###
    confs = vecNormConf(y_std)
    topConf = y_mean + confs
    botConf = y_mean - confs
    plt.plot(X, y_mean, label="Plot")
    plt.fill_between(np.reshape(X, -1), botConf, topConf, color='blue', alpha=0.1)

    plt.xlabel("X")
    _ = plt.ylabel("y")
    plt.legend()
    plt.show()


def plotSPR2(spr: GaussianProcessRegressor):
    ### Plotting Space ###
    X = np.linspace(0, 10, num=100).reshape(-1, 1)
    y = target_generator(X, add_noise=False)

    rng = np.random.RandomState(0)
    X_train = np.linspace(3, 7, num=20).reshape(-1, 1)
    y_train = target_generator(X_train, add_noise=False)

    ### Predict Using Stochastic Process ###
    spr.fit(X_train, y_train)
    y_mean, y_std = spr.predict(X, return_std=True)

    ### Ground Truth and Observations ###
    plt.plot(X, y, label="Ground Truth", color='red')
    plt.scatter(
        x=X_train[:, 0],
        y=y_train,
        color="black",
        alpha=0.4,
        label="Observations",
    )

    ### Stochastic Process ###
    confs = vecNormConf(y_std)
    topConf = y_mean + confs
    botConf = y_mean - confs
    plt.plot(X, y_mean, label="Plot")
    plt.fill_between(np.reshape(X, -1), botConf, topConf, color='blue', alpha=0.1)

    plt.xlabel("X")
    _ = plt.ylabel("y")
    plt.legend()
    plt.show()






if __name__ == '__main__':
    m_test()







