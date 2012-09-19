from numpy.testing import assert_array_almost_equal, assert_almost_equal, \
                          assert_equal
from nose.tools import assert_true
import numpy as np

from sklearn.linear_model.coordinate_descent import ElasticNet, enet_path, \
                                        elastic_net_strong_rule_active_set
from sklearn.datasets.samples_generator import make_regression
from sklearn.linear_model.cd_fast import elastic_net_kkt_violating_features


def test_enet_basic_strong_rule_filtering():
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)
    alpha = 110.0
    rho = 0.9

    basic_strong_set = elastic_net_strong_rule_active_set(X, y,
                alpha=alpha, rho=rho)

    ref_strong_set = np.array([False for x in range(X.shape[1])], dtype=bool)
    ref_strong_set[[4, 33, 35]] = True

    assert_equal(basic_strong_set, ref_strong_set)


def test_enet_sequential_strong_rule_filtering():
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)
    alphas = [110.0, 75]
    rho = 0.9

    clf = ElasticNet(alpha=alphas[0], rho=rho, precompute=False)
    clf.fit(X, y)

    sequential_strong_set = elastic_net_strong_rule_active_set(X, y,
                                alpha=alphas[1], rho=rho,
                                alpha_init=alphas[0], coef_init=clf.coef_)

    ref_sequential_strong_set = np.array([False for x in range(X.shape[1])],
                                                                 dtype=bool)
    ref_sequential_strong_set[[4, 8, 12, 14, 16, 22, 25, 29, 31, \
                                        33, 35, 37, 40, 49]] = True

    assert_equal(sequential_strong_set, ref_sequential_strong_set)


def test_automatic_strong_rule_selection():
    """
    test that sequential strong rule is used, if the needed
    information are provided
    """
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)
    alphas = [110.0, 75.0]
    rho = 0.9

    clf = ElasticNet(alpha=alphas[0], rho=rho, precompute=False)
    clf.fit(X, y)

    basic_strong_set = elastic_net_strong_rule_active_set(X, y, \
                alpha=alphas[1], rho=rho)

    sequential_strong_set = elastic_net_strong_rule_active_set(X, y, \
                                alpha=alphas[1], rho=rho, \
                                alpha_init=alphas[0], coef_init=clf.coef_)

    # the sequential strong rule is expected to select less coefs
    assert_true(np.count_nonzero(sequential_strong_set) <
                np.count_nonzero(basic_strong_set))


def test_elastic_net_find_kkt_violators():
    # watch out!!! this test only passed after fiddling with the tol of the
    # kkt check
    MAX_ITER = 1000
    X, y = make_regression(n_samples=40, n_features=10, n_informative=5,
                     random_state=0)
    alpha = 50.0
    rho = 0.9
    tol = 1e-6

    clf = ElasticNet(alpha=alpha, rho=rho, tol=tol, precompute=False)
    clf.fit(X, y)
    changed_coefs = clf.coef_

    subset = [0, 1, 2, 3, 4, 7]
    n_samples = X.shape[0]
    l1_reg = alpha * rho * n_samples
    l2_reg = alpha * (1.0 - rho) * n_samples

    # change some coefs
    changed_coefs[3] += 2
    changed_coefs[5] += 2
    changed_coefs[7] += -2

    R = y - np.dot(X, changed_coefs)
    X = np.asfortranarray(X)

    found_violators_in_subset = elastic_net_kkt_violating_features(
            changed_coefs, l1_reg, l2_reg, X, y, R, subset=subset)
    print found_violators_in_subset
    assert_true(found_violators_in_subset == set([3, 7]))

    found_violators_in_full_set = elastic_net_kkt_violating_features(
            changed_coefs, l1_reg, l2_reg, X, y, R, subset=None)
    assert_true(found_violators_in_full_set == set([3, 5, 7]))


def test_enet_strong_rule_against_standart_enet():
    X, y = make_regression(n_samples=40, n_features=20, n_informative=5,
                    random_state=0)
    alpha = 10
    rho = 0.9

    clf = ElasticNet(alpha=alpha, rho=rho)
    clf.fit(X, y)

    clf_strong_rule = ElasticNet(alpha=alpha, rho=rho, precompute=False,
                                  use_strong_rule=True)
    clf_strong_rule.fit(X, y)

    assert_array_almost_equal(clf.coef_, clf_strong_rule.coef_, 3)

    alpha = 0.5
    MAX_ITER = 1000
    tol = 1e-5

    clf = ElasticNet(alpha=alpha, rho=rho, tol=tol, max_iter=MAX_ITER)
    clf.fit(X, y)

    clf_strong_rule = ElasticNet(alpha=alpha, rho=rho, precompute=False,
                        use_strong_rule=True, tol=tol, max_iter=MAX_ITER)
    clf_strong_rule.fit(X, y)

    assert_array_almost_equal(clf.coef_, clf_strong_rule.coef_, 3)


def test_enet_path():
    MAX_ITER = 1000
    X, y = make_regression(n_samples=100, n_features=50, n_informative=20,
                     random_state=0)

    rho = 0.9
    alphas = np.array([0.2, 38, 75, 110])

    models = enet_path(X, y, tol=1e-8, alphas=alphas, rho=rho,
        fit_intercept=False, normalize=False, precompute=False, copy_X=True,
        max_iter=MAX_ITER, use_strong_rule=True)
    coefs_sr_ = np.array([m.coef_ for m in models])

    # compare with regular enet_path
    models = enet_path(X, y, tol=1e-8, alphas=alphas, rho=rho,
                       fit_intercept=False, normalize=False, precompute=False,
                       copy_X=True, max_iter=MAX_ITER, use_strong_rule=False)
    coefs_skl_ = np.array([m.coef_ for m in models])

    assert_array_almost_equal(coefs_sr_.flatten(), coefs_skl_.flatten(), 5)

    # test with automatic sequence of alphas
    MAX_ITER = 10000
    rho = 0.9

    models = enet_path(X, y, tol=1e-14, rho=rho,
        fit_intercept=False, normalize=False, precompute=False, copy_X=True,
        max_iter=MAX_ITER, use_strong_rule=True)
    coefs_sr_ = np.array([m.coef_ for m in models])

    # compare with regular enet_path
    models = enet_path(X, y, tol=1e-14, rho=rho,
                       fit_intercept=False, normalize=False, precompute=False,
                       copy_X=True, max_iter=MAX_ITER, use_strong_rule=False)
    coefs_skl_ = np.array([m.coef_ for m in models])

    assert_array_almost_equal(coefs_sr_.flatten(), coefs_skl_.flatten())
