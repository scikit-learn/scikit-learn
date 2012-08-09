from numpy.testing import assert_array_almost_equal, assert_almost_equal, assert_equal
import numpy as np

from sklearn.linear_model.coordinate_descent import ElasticNet, enet_path
from sklearn.datasets.samples_generator import make_regression


def test_enet_basic_strong_rule_filtering():
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)
    alpha = 110.0
    rho = 0.9

    clf_strong_rule = ElasticNet(alpha=alpha, rho=rho, max_iter=MAX_ITER, \
                                  precompute=False, use_strong_rule=True)

    basic_strong_set = np.array(clf_strong_rule._filter_with_strong_rule(X, y))
    assert_equal(basic_strong_set, [4, 33, 35] )


def test_enet_sequential_strong_rule_filtering():
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)
    alpha_0 = 110.0
    alpha_1 = 75.0
    rho = 0.9

    clf = ElasticNet(alpha=alpha_0, rho=rho, precompute=False)
    clf.fit(X, y)

    clf_strong_rule = ElasticNet(alpha=alpha_1, rho=rho, precompute=False, \
                                  use_strong_rule=True)
    # the sequential strong rule needs alpha_{n-1} and the corresponding coefs
    #clf_strong_rule.coef_ = clf.coef_

    sequential_strong_set = \
        np.array(clf_strong_rule._filter_with_strong_rule(X, y, \
                                 last_alpha=alpha_0, last_coef=clf.coef_))
    assert_equal(sequential_strong_set, \
                  [4, 8, 12, 14, 16, 22, 25, 29, 31, 33, 35, 37, 40, 49] )


def test_automatic_strong_rule_selection():
    """
    test that sequential strong rule is selected, if bigger alpha including
    coefs is given
    """
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)
    alphas = [110.0, 75.0]
    rho = 0.9

    clf_strong_rule = ElasticNet(alpha=alphas[1], rho=rho, \
                    max_iter=MAX_ITER, precompute=False, use_strong_rule=True)
    basic_strong_set = \
        np.array(clf_strong_rule._filter_with_strong_rule(X, y, \
                                             last_alpha=None, last_coef=None))

    clf = ElasticNet(alpha=alphas[0], rho=rho, precompute=False)
    clf.fit(X,y)

    clf_strong_rule = ElasticNet(alpha=alphas[1], rho=rho, precompute=False, use_strong_rule=True)
    sequential_strong_set = np.array(clf_strong_rule._filter_with_strong_rule(X, y, last_alpha=alphas[0], \
                                                                              last_coef=clf.coef_))
    # the sequential strong rule is expected to select less coefs
    assert(len(sequential_strong_set) < len(basic_strong_set))


def test_enet_kkt_check_on_subset():
    # watch out!!! this test only passed after fiddling with the tol of the
    # kkt check
    MAX_ITER = 1000
    X, y = make_regression(n_samples=40, n_features=10, n_informative=5,
                     random_state=0)
    alpha = 50
    rho = 0.9
    tol = 1e-6

    clf = ElasticNet(alpha=alpha, rho=rho, tol=tol, precompute=False)
    clf.fit(X, y)

    subset = [0, 1, 2, 3, 4, 7]
    n_samples = X.shape[0]
    l1_reg = alpha * rho * n_samples
    l2_reg = alpha * (1.0 - rho) * n_samples

    clf_strong_rule = ElasticNet(alpha=alpha, rho=rho, precompute=False, \
                                  use_strong_rule=True)
    clf_strong_rule.coef_ = clf.coef_.copy()

    # change two coefs
    clf_strong_rule.coef_[3] += 2
    clf_strong_rule.coef_[7] += -2
    clf_strong_rule._ever_active_set = []
    clf_strong_rule._enet_add_kkt_violating_features(X, y, l1_reg, l2_reg, \
                                     subset=subset)

    ever_active_set = np.array(clf_strong_rule._ever_active_set)
    assert_equal(ever_active_set, [3, 7])


def test_enet_kkt_check_on_full_set():
    # watch out!!! this test only passed after fiddling with the tol of the
    # kkt check
    MAX_ITER = 1000
    X, y = make_regression(n_samples=40, n_features=10, n_informative=5,
                     random_state=0)
    alpha = 50
    rho = 0.9
    tol = 1e-6

    clf = ElasticNet(alpha=alpha, rho=rho, tol=tol, precompute=False)
    clf.fit(X, y)

    n_samples = X.shape[0]
    l1_reg = alpha * rho * n_samples
    l2_reg = alpha * (1.0 - rho) * n_samples

    clf_strong_rule = ElasticNet(alpha=alpha, rho=rho, precompute=False, \
                                  use_strong_rule=True)
    clf_strong_rule.coef_ = clf.coef_.copy()

    # change two coefs
    clf_strong_rule.coef_[3] += 2
    clf_strong_rule.coef_[7] += -2
    clf_strong_rule._ever_active_set = []
    clf_strong_rule._enet_add_kkt_violating_features(X, y, l1_reg, l2_reg, \
                                     subset=None)

    ever_active_set = np.array(clf_strong_rule._ever_active_set)
    assert_equal(ever_active_set, [3, 7])


def test_enet_strong_rule_against_standart_enet():
    X, y = make_regression(n_samples=40, n_features=20, n_informative=5,
                    random_state=0)
    alpha = 10
    rho = 0.9

    clf = ElasticNet(alpha=alpha, rho=rho)
    clf.fit(X, y)

    clf_strong_rule = ElasticNet(alpha=alpha, rho=rho, precompute=False, \
                                  use_strong_rule=True)
    clf_strong_rule.fit(X, y)

    assert_array_almost_equal(clf.coef_, clf_strong_rule.coef_)
#    assert_equal(0,1)


def test_enet_path():
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)

    rho = 0.9
    alphas = np.array([0.2, 38, 75, 110])

    models = enet_path(X, y, tol=1e-8, alphas=alphas, rho=rho, \
        fit_intercept=False, normalize=False, precompute=False, copy_X=True, \
             max_iter=MAX_ITER, use_strong_rule=True)
    coefs_sr_ = np.array([m.coef_ for m in models])

    # compare with regular enet_path
    models = enet_path(X, y, tol=1e-8, alphas=alphas, rho=rho, \
                       fit_intercept=False, normalize=False, precompute=False, \
                       copy_X=True, max_iter=MAX_ITER, use_strong_rule=False)
    coefs_skl_ = np.array([m.coef_ for m in models])

    # flatten array and compare
    assert_array_almost_equal(coefs_sr_[1], coefs_skl_[1], 5)
    assert_array_almost_equal(coefs_sr_[3], coefs_skl_[3], 5)
#    assert_equal(0,1)