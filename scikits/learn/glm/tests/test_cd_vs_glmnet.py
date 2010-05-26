# Authors: Olivier Grisel <olivier.grisel@ensta.org>
#          Alexandre Gramfort <alexandre.gramfort@inria.fr>
# License: BSD Style.

# $Id$

import numpy as np
import scipy.linalg as linalg

from numpy.testing import assert_array_almost_equal

from ..coordinate_descent import Lasso, lasso_path
from ..coordinate_descent import ElasticNet, enet_path

from glmnet.elastic_net import Lasso as LassoGLMNET
from glmnet.elastic_net import ElasticNet as ElasticNetGLMNET
from parietal.learn.sparse.elastic_net import Lasso as LassoProx
from parietal.learn.sparse.elastic_net import ElasticNet as ElasticNetProx

def lasso_objective(X, y, w, alpha):
    R = y - np.dot(X, w)
    return 0.5 / X.shape[0] * linalg.norm(R) ** 2 + alpha * np.abs(w).sum()

def enet_objective(X, y, w, alpha, rho):
    R = y - np.dot(X, w)
    return 0.5 / X.shape[0] * linalg.norm(R) ** 2 + rho * alpha * np.abs(w).sum() \
            + 0.5 * alpha * (1-rho) * (w**2).sum()

def test_lasso_cd_vs_glmnet():
    """
    test scikit lasso based on coordinate descent agains the results
    obtained with GLMNET

    """

    # generate toy dataset
    np.random.seed(0)
    X = np.random.randn(3, 5)
    w = 3*np.random.randn(5)
    w[2:] = 0
    y = np.dot(X, w)

    # remove intercept by standardizing data
    y -= np.mean(y)
    X -= X.mean(0)
    X /= np.sqrt(np.sum(X**2, 0))

    alpha = 0.1

    # CD
    lasso = Lasso(alpha=alpha, tol=1e-18)
    y_pred = lasso.fit(X, y, maxit=10000).predict(X)

    # GLMNET
    lasso_glmnet = LassoGLMNET(alpha=alpha)
    y_pred_glmnet = lasso_glmnet.fit(X, y, threshold=1e-15).predict(X)

    # PROX
    lasso_prox = LassoProx(alpha=alpha)
    y_pred_prox = lasso_prox.fit(X, y, maxit=5000, tol=1e-12).predict(X)

    print "-- objective"
    print lasso_objective(X, y, lasso.coef_, alpha), ' CD'
    print lasso_objective(X, y, lasso_glmnet.coef_.ravel(), alpha), ' GLMNET'
    print lasso_objective(X, y, lasso_prox.coef_.ravel(), alpha), ' PROX'

    print "-- KKT"
    # print alpha * X.shape[0]
    print np.dot(X.T, y - np.dot(X, lasso.coef_)), ' CD'
    print np.dot(X.T, y - np.dot(X, lasso_glmnet.coef_.ravel())), ' GLMNET'
    print np.dot(X.T, y - np.dot(X, lasso_prox.coef_.ravel())), ' PROX'

    print "-- coefs"
    print w , ' TRUE'
    print lasso.coef_, ' CD'
    print lasso_glmnet.coef_.ravel().T, ' GLMNET'
    print lasso_prox.coef_.ravel().T, ' PROX'
    print "-- intercept"
    print lasso_glmnet.intercept_

def test_enet_cd_vs_glmnet():
    """
    test scikit lasso based on coordinate descent agains the results
    obtained with GLMNET

    """

    # generate toy dataset
    np.random.seed(0)
    X = np.random.randn(3, 5)
    w = 3*np.random.randn(5)
    w[2:] = 0
    y = np.dot(X, w)

    # remove intercept by standardizing data
    y -= np.mean(y)
    X -= X.mean(0)
    X /= np.sqrt(np.sum(X**2, 0))

    alpha = 1
    rho = 0.5

    # CD
    enet = ElasticNet(alpha=alpha, rho=rho, tol=1e-12)
    y_pred = enet.fit(X, y, maxit=10000).predict(X)

    # GLMNET
    enet_glmnet = ElasticNetGLMNET(alpha=alpha, rho=rho)
    y_pred_glmnet = enet_glmnet.fit(X, y, threshold=1e-15).predict(X)

    # PROX
    enet_prox = ElasticNetProx(alpha=alpha, rho=rho)
    y_pred_prox = enet_prox.fit(X, y, maxit=5000, tol=1e-12).predict(X)

    print "-- objective"
    print enet_objective(X, y, enet.coef_, alpha, rho), ' CD'
    print enet_objective(X, y, enet_glmnet.coef_.ravel(), alpha, rho), ' GLMNET'
    print enet_objective(X, y, enet_prox.coef_.ravel(), alpha, rho), ' PROX'

    print "-- KKT"
    # print alpha * X.shape[0]
    print np.dot(X.T, y - np.dot(X, enet.coef_)), ' CD'
    print np.dot(X.T, y - np.dot(X, enet_glmnet.coef_.ravel())), ' GLMNET'
    print np.dot(X.T, y - np.dot(X, enet_prox.coef_.ravel())), ' PROX'

    print "-- coefs"
    print w , ' TRUE'
    print enet.coef_, ' CD'
    print enet_glmnet.coef_.ravel().T, ' GLMNET'
    print enet_prox.coef_.ravel().T, ' PROX'
    print "-- intercept"
    print enet_glmnet.intercept_

