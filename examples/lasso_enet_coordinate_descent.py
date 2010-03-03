# Author: Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

from itertools import cycle
import numpy as np
import pylab as pl
from scikits.learn.glm.cd import Lasso, ElasticNet, lasso_path, enet_path

def demo_glm():
    """
    Sample usage of GLM to predict with a linear model.
    It will plot the paths for the Lasso and the Elastic-Net problem

    """
    n_samples, n_features, maxit = 5, 10, 30
    np.random.seed(0)
    y = np.random.randn(n_samples)
    X = np.random.randn(n_samples, n_features)

    """Lasso model
    """
    alpha = 1.0

    lasso = Lasso(alpha=alpha)
    lasso.fit(X, y, maxit=maxit)

    print "Duality gap Lasso (should be small): %f" % lasso.gap

    import pylab as pl
    pl.close('all')
    pl.plot(lasso.E)
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.title('Lasso')
    pl.show()

    """Elastic-Net model
    """
    alpha = 1.0
    beta = 1.0

    enet = ElasticNet(alpha=alpha, beta=beta)
    enet.fit(X, y, maxit=maxit)

    print "Duality gap (should be small): %f" % enet.gap

    pl.figure()
    pl.plot(enet.E)
    pl.xlabel('Iteration')
    pl.ylabel('Cost function')
    pl.title('Elastic-Net')
    pl.show()

    """Test path functions
    """

    alphas_lasso, weights_lasso = lasso_path(X, y, factor=0.97, n_alphas = 100)
    alphas_enet, weights_enet = enet_path(X, y, factor=0.97, n_alphas = 100, beta=0.1)

    pl.figure()
    color_iter = cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    for color, weight_lasso, weight_enet in zip(color_iter, weights_lasso.T, weights_enet.T):
        pl.plot(-np.log(alphas_lasso), weight_lasso, color)
        pl.plot(-np.log(alphas_enet), weight_enet, color+'x')
    pl.xlabel('-log(lambda)')
    pl.ylabel('weights')
    pl.title('Lasso and Elastic-Net Paths')
    pl.legend(['Lasso','Elastic-Net'])
    pl.show()

if __name__ == '__main__':
    demo_glm()
