"""Testing for T process regression"""

# Authors: Conrad Stevens <conrad.stevens@sydney.edu.au>
# Modified by: TODO
# License: TODO

import re
import sys
import warnings

import numpy as np
import pytest
from scipy.optimize import approx_fprime

from sklearn.exceptions import ConvergenceWarning
from sklearn.gaussian_process import TProcessRegressor
from sklearn.gaussian_process.kernels import (
    RBF,
    DotProduct,
    ExpSineSquared,
    WhiteKernel,
)
from sklearn.gaussian_process.kernels import (
    ConstantKernel as C,
)
from sklearn.gaussian_process.tests._mini_sequence_kernel import MiniSeqKernel
from sklearn.utils._testing import (
    assert_allclose,
    assert_almost_equal,
    assert_array_almost_equal,
    assert_array_less,
)


def f(x):
    return x * np.sin(x)


X = np.atleast_2d([1.0, 3.0, 5.0, 6.0, 7.0, 8.0]).T
X2 = np.atleast_2d([2.0, 4.0, 5.5, 6.5, 7.5]).T
y = f(X).ravel()

fixed_kernel = RBF(length_scale=1.0, length_scale_bounds="fixed")
kernels = [
    RBF(length_scale=1.0),
    fixed_kernel,
    RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3)),
    C(1.0, (1e-2, 1e2)) * RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3)),
    C(1.0, (1e-2, 1e2)) * RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3))
    + C(1e-5, (1e-5, 1e2)),
    C(0.1, (1e-2, 1e2)) * RBF(length_scale=1.0, length_scale_bounds=(1e-3, 1e3))
    + C(1e-5, (1e-5, 1e2)),
]
non_fixed_kernels = [kernel for kernel in kernels if kernel != fixed_kernel]


@pytest.mark.parametrize("kernel", kernels)
def test_tpr_prediction(kernel):
    if sys.maxsize <= 2**32:
        pytest.xfail("This test may fail on 32 bit Python")

    # Test the interpolating property for different kernels.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    y_pred, y_cov = tpr.predict(X, return_cov=True)

    assert_almost_equal(y_pred, y)
    assert_almost_equal(np.diag(y_cov), 0.0)


@pytest.mark.parametrize("kernel", kernels)
def test_tpr_interpolation(kernel):
    if sys.maxsize <= 2**32:
        pytest.xfail("This test may fail on 32 bit Python")

    # Test the interpolating property for different kernels.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    y_pred, y_cov = tpr.predict(X, return_cov=True)

    assert_almost_equal(y_pred, y)
    assert_almost_equal(np.diag(y_cov), 0.0)


def test_tpr_interpolation_structured():
    # Test the interpolating property for different kernels.
    kernel = MiniSeqKernel(baseline_similarity_bounds="fixed")
    X = ["A", "B", "C"]
    y = np.array([1, 2, 3])
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    y_pred, y_cov = tpr.predict(X, return_cov=True)

    assert_almost_equal(
        kernel(X, eval_gradient=True)[1].ravel(), (1 - np.eye(len(X))).ravel()
    )
    assert_almost_equal(y_pred, y)
    assert_almost_equal(np.diag(y_cov), 0.0)


@pytest.mark.parametrize("kernel", non_fixed_kernels)
def test_lml_improving(kernel):
    if sys.maxsize <= 2**32:
        pytest.xfail("This test may fail on 32 bit Python")

    # Test that hyperparameter-tuning improves log-marginal likelihood.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    assert tpr.log_marginal_likelihood(tpr.kernel_.theta) > tpr.log_marginal_likelihood(
        kernel.theta
    )


@pytest.mark.parametrize("kernel", kernels)
def test_lml_precomputed(kernel):
    # Test that lml of optimized kernel is stored correctly.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    assert tpr.log_marginal_likelihood(tpr.kernel_.theta) == pytest.approx(
        tpr.log_marginal_likelihood()
    )


@pytest.mark.parametrize("kernel", kernels)
def test_lml_without_cloning_kernel(kernel):
    # Test that lml of optimized kernel is stored correctly.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    input_theta = np.ones(tpr.kernel_.theta.shape, dtype=np.float64)

    tpr.log_marginal_likelihood(input_theta, clone_kernel=False)
    assert_almost_equal(tpr.kernel_.theta, input_theta, 7)


@pytest.mark.parametrize("kernel", non_fixed_kernels)
def test_converged_to_local_maximum(kernel):
    # Test that we are in local maximum after hyperparameter-optimization.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)

    lml, lml_gradient = tpr.log_marginal_likelihood(tpr.kernel_.theta, True)

    assert np.all(
        (np.abs(lml_gradient) < 1e-4)
        | (tpr.kernel_.theta == tpr.kernel_.bounds[:, 0])
        | (tpr.kernel_.theta == tpr.kernel_.bounds[:, 1])
    )


@pytest.mark.parametrize("kernel", non_fixed_kernels)
def test_solution_inside_bounds(kernel):
    # Test that hyperparameter-optimization remains in bounds#
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)

    bounds = tpr.kernel_.bounds
    max_ = np.finfo(tpr.kernel_.theta.dtype).max
    tiny = 1e-10
    bounds[~np.isfinite(bounds[:, 1]), 1] = max_

    assert_array_less(bounds[:, 0], tpr.kernel_.theta + tiny)
    assert_array_less(tpr.kernel_.theta, bounds[:, 1] + tiny)


@pytest.mark.parametrize("kernel", kernels)
def test_lml_gradient(kernel):
    # Compare analytic and numeric gradient of log marginal likelihood.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)

    lml, lml_gradient = tpr.log_marginal_likelihood(kernel.theta, True)
    lml_gradient_approx = approx_fprime(
        kernel.theta, lambda theta: tpr.log_marginal_likelihood(theta, False), 1e-10
    )

    assert_almost_equal(lml_gradient, lml_gradient_approx, 3)


@pytest.mark.parametrize("kernel", kernels)
def test_prior(kernel):
    # Test that TP prior has mean 0 and identical variances.
    tpr = TProcessRegressor(kernel=kernel)

    y_mean, y_cov = tpr.predict(X, return_cov=True)

    assert_almost_equal(y_mean, 0, 5)
    if len(tpr.kernel.theta) > 1:
        # XXX: quite hacky, works only for current kernels
        assert_almost_equal(np.diag(y_cov), np.exp(kernel.theta[0]), 5)
    else:
        assert_almost_equal(np.diag(y_cov), 1, 5)


@pytest.mark.parametrize("kernel", kernels)
def test_sample_statistics(kernel):
    # Test that statistics of samples drawn from TP are correct.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)

    y_mean, y_cov = tpr.predict(X2, return_cov=True)

    samples = tpr.sample_y(X2, 300000)

    # More digits accuracy would require many more samples
    assert_almost_equal(y_mean, np.mean(samples, 1), 1)
    assert_almost_equal(
        np.diag(y_cov) / np.diag(y_cov).max(),
        np.var(samples, 1) / np.diag(y_cov).max(),
        1,
    )


def test_no_optimizer():
    # Test that kernel parameters are unmodified when optimizer is None.
    kernel = RBF(1.0)
    tpr = TProcessRegressor(kernel=kernel, optimizer=None).fit(X, y)
    assert np.exp(tpr.kernel_.theta) == 1.0


@pytest.mark.parametrize("kernel", kernels)
@pytest.mark.parametrize("target", [y, np.ones(X.shape[0], dtype=np.float64)])
def test_predict_cov_vs_std(kernel, target):
    if sys.maxsize <= 2**32:
        pytest.xfail("This test may fail on 32 bit Python")

    # Test that predicted std.-dev. is consistent with cov's diagonal.
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    y_mean, y_cov = tpr.predict(X2, return_cov=True)
    y_mean, y_std = tpr.predict(X2, return_std=True)
    assert_almost_equal(np.sqrt(np.diag(y_cov)), y_std)


def test_anisotropic_kernel():
    # Test that tpr can identify meaningful anisotropic length-scales.
    # We learn a function which varies in one dimension ten-times slower
    # than in the other. The corresponding length-scales should differ by at
    # least a factor 5
    rng = np.random.RandomState(0)
    # X = rng.uniform(-1, 1, (50, 2))
    X = rng.uniform(
        -1, 1, (25, 2)
    )  # Less observations must be used for testing T-process Regression
    y = X[:, 0] + 0.1 * X[:, 1]

    kernel = RBF([1.0, 1.0])
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    assert np.exp(tpr.kernel_.theta[1]) > np.exp(tpr.kernel_.theta[0]) * 5


def test_random_starts():
    # Test that an increasing number of random-starts of TP fitting only
    # increases the log marginal likelihood of the chosen theta.
    n_samples, n_features = 25, 2
    rng = np.random.RandomState(0)
    X = rng.randn(n_samples, n_features) * 2 - 1
    y = (
        np.sin(X).sum(axis=1)
        + np.sin(3 * X).sum(axis=1)
        + rng.normal(scale=0.1, size=n_samples)
    )

    kernel = C(1.0, (1e-2, 1e2)) * RBF(
        length_scale=[1.0] * n_features, length_scale_bounds=[(1e-4, 1e2)] * n_features
    ) + WhiteKernel(noise_level=1e-5, noise_level_bounds=(1e-5, 1e1))
    last_lml = -np.inf
    # for n_restarts_optimizer in range(5):
    for n_restarts_optimizer in range(
        3
    ):  # Less restarts must be used when testing T-Processes
        tp = TProcessRegressor(
            kernel=kernel,
            n_restarts_optimizer=n_restarts_optimizer,
            random_state=0,
        ).fit(X, y)
        lml = tp.log_marginal_likelihood(tp.kernel_.theta)
        assert lml >= last_lml - np.finfo(np.float32).eps  # Made >= in the TP case
        last_lml = lml


@pytest.mark.parametrize("kernel", kernels)
def test_y_normalization(kernel):
    """
    Test normalization of the target values in TP

    Fitting non-normalizing TP on normalized y and fitting normalizing TP
    on unnormalized y should yield identical results. Note that, here,
    'normalized y' refers to y that has been made zero mean and unit
    variance.

    """

    y_mean = np.mean(y)
    y_std = np.std(y)
    y_norm = (y - y_mean) / y_std

    # Fit non-normalizing TP on normalized y
    tpr = TProcessRegressor(kernel=kernel)
    tpr.fit(X, y_norm)

    # Fit normalizing TP on unnormalized y
    tpr_norm = TProcessRegressor(kernel=kernel, normalize_y=True)
    tpr_norm.fit(X, y)

    # Compare predicted mean, std-devs and covariances
    y_pred, y_pred_std = tpr.predict(X2, return_std=True)
    y_pred = y_pred * y_std + y_mean
    y_pred_std = y_pred_std * y_std
    y_pred_norm, y_pred_std_norm = tpr_norm.predict(X2, return_std=True)

    assert_almost_equal(y_pred, y_pred_norm)
    assert_almost_equal(y_pred_std, y_pred_std_norm)

    _, y_cov = tpr.predict(X2, return_cov=True)
    y_cov = y_cov * y_std**2
    _, y_cov_norm = tpr_norm.predict(X2, return_cov=True)

    assert_almost_equal(y_cov, y_cov_norm)


def test_large_variance_y():
    """
    Here we test that, when noramlize_y=True, our TP (acting as a GP)
    can produce a sensible fit to training data whose variance is
    significantly larger than unity. This test was made in response
    to issue #15612.

    TP predictions with large starting degrees of freedom, making them
    act as a GP are verified against predictions that were made
    using GPy which, here, is treated as the 'gold standard'. Note that we
    only investigate the RBF kernel here, as that is what was used in the
    GPy implementation.

    Note, the mean and standard deviation should still be the same as
    given by a gp because the kernel parameters are being trained on.
    The marginal distributions at each point are still TPs though.

    The following code can be used to recreate the GPy data:

    --------------------------------------------------------------------------
    import GPy

    kernel_gpy = GPy.kern.RBF(input_dim=1, lengthscale=1.)
    gpy = GPy.models.tpregression(X, np.vstack(y_large), kernel_gpy)
    gpy.optimize()
    y_pred_gpy, y_var_gpy = gpy.predict(X2)
    y_pred_std_gpy = np.sqrt(y_var_gpy)
    --------------------------------------------------------------------------
    """

    # Here we utilise a larger variance version of the training data
    y_large = 10 * y

    # Standard TP with normalize_y=True
    RBF_params = {"length_scale": 1.0}
    kernel = RBF(**RBF_params)
    tpr = TProcessRegressor(kernel=kernel, normalize_y=True, v=1_000)
    tpr.fit(X, y_large)
    y_pred, y_pred_std = tpr.predict(X2, return_std=True)

    # 'Gold standard' mean predictions from GPy
    y_pred_gpy = np.array(
        [15.16918303, -27.98707845, -39.31636019, 14.52605515, 69.18503589]
    )

    # 'Gold standard' std predictions from GPy
    y_pred_std_gpy = np.array(
        [7.78860962, 3.83179178, 0.63149951, 0.52745188, 0.86170042]
    )

    # Based on numerical experiments, it's reasonable to expect our
    # TP's mean predictions to get within 7% of predictions of those
    # made by GPy.
    assert_allclose(y_pred, y_pred_gpy, rtol=0.07, atol=0)

    # Based on numerical experiments, it's reasonable to expect our
    # TP's std predictions to get within 15% of predictions of those
    # made by GPy.
    assert_allclose(y_pred_std, y_pred_std_gpy, rtol=0.15, atol=0)


def test_y_multioutput():
    # Test that tpr_1 can deal with multi-dimensional target values
    y_2d = np.vstack((y, y * 2)).T

    # Test for fixed kernel that first dimension of 2d TP equals the output
    # of 1d TP and that second dimension is twice as large
    kernel = RBF(length_scale=1.0)

    tpr_1 = TProcessRegressor(kernel=kernel, optimizer=None, normalize_y=False)
    tpr_2 = TProcessRegressor(kernel=kernel, optimizer=None, normalize_y=False)
    tpr_1.fit(X, y)
    tpr_2.fit(X, y * 2)

    tpr_2d = TProcessRegressor(kernel=kernel, optimizer=None, normalize_y=False)
    tpr_2d.fit(X, y_2d)

    y_pred_1d, y_1_std_1d = tpr_1.predict(X2, return_std=True)
    _, y_2_std_1d = tpr_2.predict(X2, return_std=True)
    y_pred_2d, y_std_2d = tpr_2d.predict(X2, return_std=True)

    _, y_1_cov_1d = tpr_1.predict(X2, return_cov=True)
    _, y_2_cov_1d = tpr_2.predict(X2, return_cov=True)
    _, y_cov_2d = tpr_2d.predict(X2, return_cov=True)

    assert_almost_equal(y_pred_1d, y_pred_2d[:, 0])
    assert_almost_equal(y_pred_1d, y_pred_2d[:, 1] / 2)

    # Standard deviation and covariance do not depend on output
    assert_almost_equal(y_1_std_1d, y_std_2d[..., 0])
    assert_almost_equal(y_2_std_1d, y_std_2d[..., 1])

    y_sample_1d = tpr_1.sample_y(X2, n_samples=10)
    y_sample_2d = tpr_2d.sample_y(X2, n_samples=10)

    assert y_sample_1d.shape == (5, 10)
    assert y_sample_2d.shape == (5, 2, 10)
    # Only the first target will be equal
    assert_almost_equal(y_sample_1d, y_sample_2d[:, 0, :])

    # Test hyperparameter optimization
    for kernel in kernels:
        tpr_1 = TProcessRegressor(kernel=kernel, normalize_y=True)
        tpr_1.fit(X, y)

        tpr_2d = TProcessRegressor(kernel=kernel, normalize_y=True)
        tpr_2d.fit(X, np.vstack((y, y)).T)

        assert_almost_equal(tpr_1.kernel_.theta, tpr_2d.kernel_.theta, 4)


@pytest.mark.parametrize("kernel", non_fixed_kernels)
def test_custom_optimizer(kernel):
    # Test that tpr can use externally defined optimizers.
    # Define a dummy optimizer that simply tests 50 random hyperparameters
    def optimizer(obj_func, initial_theta, bounds):
        rng = np.random.RandomState(0)
        theta_opt, func_min = initial_theta, obj_func(
            initial_theta, eval_gradient=False
        )
        for _ in range(50):
            theta = np.atleast_1d(
                rng.uniform(np.maximum(-2, bounds[:, 0]), np.minimum(1, bounds[:, 1]))
            )
            f = obj_func(theta, eval_gradient=False)
            if f < func_min:
                theta_opt, func_min = theta, f
        return theta_opt, func_min

    tpr = TProcessRegressor(kernel=kernel, optimizer=optimizer)
    tpr.fit(X, y)
    # Checks that optimizer improved marginal likelihood
    assert tpr.log_marginal_likelihood(tpr.kernel_.theta) > tpr.log_marginal_likelihood(
        tpr.kernel.theta
    )


def test_tpr_correct_error_message():
    X = np.arange(12).reshape(6, -1)
    y = np.ones(6)
    kernel = DotProduct()
    tpr = TProcessRegressor(kernel=kernel, alpha=0.0)
    message = (
        "The kernel, %s, is not returning a "
        "positive definite matrix. Try gradually increasing "
        "the 'alpha' parameter of your "
        "GaussianProcessRegressor estimator." % kernel
    )
    with pytest.raises(np.linalg.LinAlgError, match=re.escape(message)):
        tpr.fit(X, y)


@pytest.mark.parametrize("kernel", kernels)
def test_duplicate_input(kernel):
    # Test tpr can handle two different output-values for the same input.
    tpr_equal_inputs = TProcessRegressor(kernel=kernel, alpha=1e-2)
    tpr_similar_inputs = TProcessRegressor(kernel=kernel, alpha=1e-2)

    X_ = np.vstack((X, X[0]))
    y_ = np.hstack((y, y[0] + 1))
    tpr_equal_inputs.fit(X_, y_)

    X_ = np.vstack((X, X[0] + 1e-15))
    y_ = np.hstack((y, y[0] + 1))
    tpr_similar_inputs.fit(X_, y_)

    X_test = np.linspace(0, 10, 100)[:, None]
    y_pred_equal, y_std_equal = tpr_equal_inputs.predict(X_test, return_std=True)
    y_pred_similar, y_std_similar = tpr_similar_inputs.predict(X_test, return_std=True)

    assert_almost_equal(y_pred_equal, y_pred_similar)
    assert_almost_equal(y_std_equal, y_std_similar)


def test_no_fit_default_predict():
    # Test that tpr predictions without fit does not break by default.
    default_kernel = C(1.0, constant_value_bounds="fixed") * RBF(
        1.0, length_scale_bounds="fixed"
    )
    tpr1 = TProcessRegressor()
    _, y_std1 = tpr1.predict(X, return_std=True)
    _, y_cov1 = tpr1.predict(X, return_cov=True)

    tpr2 = TProcessRegressor(kernel=default_kernel)
    _, y_std2 = tpr2.predict(X, return_std=True)
    _, y_cov2 = tpr2.predict(X, return_cov=True)

    assert_array_almost_equal(y_std1, y_std2)
    assert_array_almost_equal(y_cov1, y_cov2)


def test_warning_bounds():
    kernel = RBF(length_scale_bounds=[1e-5, 1e-3])
    tpr = TProcessRegressor(kernel=kernel)
    warning_message = (
        "The optimal value found for dimension 0 of parameter "
        "length_scale is close to the specified upper bound "
        "0.001. Increasing the bound and calling fit again may "
        "find a better value."
    )
    with pytest.warns(ConvergenceWarning, match=warning_message):
        tpr.fit(X, y)

    kernel_sum = WhiteKernel(noise_level_bounds=[1e-5, 1e-3]) + RBF(
        length_scale_bounds=[1e3, 1e5]
    )
    tpr_sum = TProcessRegressor(kernel=kernel_sum)
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        tpr_sum.fit(X, y)

        assert len(record) == 2

        assert issubclass(record[0].category, ConvergenceWarning)
        assert (
            record[0].message.args[0] == "The optimal value found for "
            "dimension 0 of parameter "
            "k1__noise_level is close to the "
            "specified upper bound 0.001. "
            "Increasing the bound and calling "
            "fit again may find a better value."
        )

        assert issubclass(record[1].category, ConvergenceWarning)
        assert (
            record[1].message.args[0] == "The optimal value found for "
            "dimension 0 of parameter "
            "k2__length_scale is close to the "
            "specified lower bound 1000.0. "
            "Decreasing the bound and calling "
            "fit again may find a better value."
        )

    X_tile = np.tile(X, 2)
    kernel_dims = RBF(length_scale=[1.0, 2.0], length_scale_bounds=[1e1, 1e2])
    tpr_dims = TProcessRegressor(kernel=kernel_dims)

    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        tpr_dims.fit(X_tile, y)

        assert len(record) == 2

        assert issubclass(record[0].category, ConvergenceWarning)
        assert (
            record[0].message.args[0] == "The optimal value found for "
            "dimension 0 of parameter "
            "length_scale is close to the "
            "specified lower bound 10.0. "
            "Decreasing the bound and calling "
            "fit again may find a better value."
        )

        assert issubclass(record[1].category, ConvergenceWarning)
        assert (
            record[1].message.args[0] == "The optimal value found for "
            "dimension 1 of parameter "
            "length_scale is close to the "
            "specified lower bound 10.0. "
            "Decreasing the bound and calling "
            "fit again may find a better value."
        )


def test_bound_check_fixed_hyperparameter():
    # Regression test for issue #17943
    # Check that having a hyperparameter with fixed bounds doesn't cause an
    # error
    k1 = 50.0**2 * RBF(length_scale=50.0)  # long term smooth rising trend
    k2 = ExpSineSquared(
        length_scale=1.0, periodicity=1.0, periodicity_bounds="fixed"
    )  # seasonal component
    kernel = k1 + k2
    TProcessRegressor(kernel=kernel).fit(X, y)


@pytest.mark.parametrize("kernel", kernels)
def test_constant_target(kernel):
    """Check that the std. dev. is affected to 1 when normalizing a constant
    feature.
    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/18318
    NaN where affected to the target when scaling due to null std. dev. with
    constant target.
    """
    y_constant = np.ones(X.shape[0], dtype=np.float64)

    tpr = TProcessRegressor(kernel=kernel, normalize_y=True)
    tpr.fit(X, y_constant)
    assert tpr._y_train_std == pytest.approx(1.0)

    y_pred, y_cov = tpr.predict(X, return_cov=True)
    assert_allclose(y_pred, y_constant)
    # set atol because we compare to zero
    assert_allclose(np.diag(y_cov), 0.0, atol=1e-9)

    # Test multi-target data
    n_samples, n_targets = X.shape[0], 2
    rng = np.random.RandomState(0)
    y = np.concatenate(
        [
            rng.normal(size=(n_samples, 1)),  # non-constant target
            np.full(shape=(n_samples, 1), fill_value=2),  # constant target
        ],
        axis=1,
    )

    tpr.fit(X, y)
    Y_pred, Y_cov = tpr.predict(X, return_cov=True)

    assert_allclose(Y_pred[:, 1], 2)
    assert_allclose(np.diag(Y_cov[..., 1]), 0.0, atol=1e-9)

    assert Y_pred.shape == (n_samples, n_targets)
    assert Y_cov.shape == (n_samples, n_samples, n_targets)


def test_tpr_consistency_std_cov_non_invertible_kernel():
    """Check the consistency between the returned std. dev. and the covariance.
    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/19936
    Inconsistencies were observed when the kernel cannot be inverted (or
    numerically stable).
    """
    kernel = C(8.98576054e05, (1e-12, 1e12)) * RBF(
        [5.91326520e02, 1.32584051e03], (1e-12, 1e12)
    ) + WhiteKernel(noise_level=1e-5)
    tpr = TProcessRegressor(kernel=kernel, alpha=0, optimizer=None)
    X_train = np.array(
        [
            [0.0, 0.0],
            [1.54919334, -0.77459667],
            [-1.54919334, 0.0],
            [0.0, -1.54919334],
            [0.77459667, 0.77459667],
            [-0.77459667, 1.54919334],
        ]
    )
    y_train = np.array(
        [
            [-2.14882017e-10],
            [-4.66975823e00],
            [4.01823986e00],
            [-1.30303674e00],
            [-1.35760156e00],
            [3.31215668e00],
        ]
    )
    tpr.fit(X_train, y_train)
    X_test = np.array(
        [
            [-1.93649167, -1.93649167],
            [1.93649167, -1.93649167],
            [-1.93649167, 1.93649167],
            [1.93649167, 1.93649167],
        ]
    )
    pred1, std = tpr.predict(X_test, return_std=True)
    pred2, cov = tpr.predict(X_test, return_cov=True)
    assert_allclose(std, np.sqrt(np.diagonal(cov)), rtol=1e-5)


@pytest.mark.parametrize(
    "params, TypeError, err_msg",
    [
        (
            {"alpha": np.zeros(100)},
            ValueError,
            "alpha must be a scalar or an array with same number of entries as y",
        ),
        (
            {
                "kernel": WhiteKernel(noise_level_bounds=(-np.inf, np.inf)),
                "n_restarts_optimizer": 2,
            },
            ValueError,
            "requires that all bounds are finite",
        ),
    ],
)
def test_tpr_fit_error(params, TypeError, err_msg):
    """Check that expected error are raised during fit."""
    tpr = TProcessRegressor(**params)
    with pytest.raises(TypeError, match=err_msg):
        tpr.fit(X, y)


def test_tpr_lml_error():
    """Check that we raise the proper error in the LML method."""
    tpr = TProcessRegressor(kernel=RBF()).fit(X, y)

    err_msg = "Gradient can only be evaluated for theta!=None"
    with pytest.raises(ValueError, match=err_msg):
        tpr.log_marginal_likelihood(eval_gradient=True)


def test_tpr_predict_error():
    """Check that we raise the proper error during predict."""
    tpr = TProcessRegressor(kernel=RBF()).fit(X, y)

    err_msg = (
        "At most one of return_std, return_cov, return_tShape "
        "or return_tShapeMatrix can be requested."
    )
    with pytest.raises(RuntimeError, match=err_msg):
        tpr.predict(X, return_cov=True, return_std=True)


@pytest.mark.parametrize("normalize_y", [True, False])
@pytest.mark.parametrize("n_targets", [None, 1, 10])
def test_predict_shapes(normalize_y, n_targets):
    """Check the shapes of y_mean, y_std, and y_cov in single-output
    (n_targets=None) and multi-output settings, including the edge case when
    n_targets=1, where the sklearn convention is to squeeze the predictions.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/17394
    https://github.com/scikit-learn/scikit-learn/issues/18065
    https://github.com/scikit-learn/scikit-learn/issues/22174
    """
    rng = np.random.RandomState(1234)

    n_features, n_samples_train, n_samples_test = 6, 9, 7

    y_train_shape = (n_samples_train,)
    if n_targets is not None:
        y_train_shape = y_train_shape + (n_targets,)

    # By convention single-output data is squeezed upon prediction
    y_test_shape = (n_samples_test,)
    if n_targets is not None and n_targets > 1:
        y_test_shape = y_test_shape + (n_targets,)

    X_train = rng.randn(n_samples_train, n_features)
    X_test = rng.randn(n_samples_test, n_features)
    y_train = rng.randn(*y_train_shape)

    model = TProcessRegressor(normalize_y=normalize_y)
    model.fit(X_train, y_train)

    y_pred, y_std = model.predict(X_test, return_std=True)
    _, y_cov = model.predict(X_test, return_cov=True)

    assert y_pred.shape == y_test_shape
    assert y_std.shape == y_test_shape
    assert y_cov.shape == (n_samples_test,) + y_test_shape


@pytest.mark.parametrize("normalize_y", [True, False])
@pytest.mark.parametrize("n_targets", [None, 1, 10])
def test_sample_y_shapes(normalize_y, n_targets):
    """Check the shapes of y_samples in single-output (n_targets=0) and
    multi-output settings, including the edge case when n_targets=1, where the
    sklearn convention is to squeeze the predictions.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/22175
    """
    rng = np.random.RandomState(1234)

    n_features, n_samples_train = 6, 9
    # Number of spatial locations to predict at
    n_samples_X_test = 7
    # Number of sample predictions per test point
    n_samples_y_test = 5

    y_train_shape = (n_samples_train,)
    if n_targets is not None:
        y_train_shape = y_train_shape + (n_targets,)

    # By convention single-output data is squeezed upon prediction
    if n_targets is not None and n_targets > 1:
        y_test_shape = (n_samples_X_test, n_targets, n_samples_y_test)
    else:
        y_test_shape = (n_samples_X_test, n_samples_y_test)

    X_train = rng.randn(n_samples_train, n_features)
    X_test = rng.randn(n_samples_X_test, n_features)
    y_train = rng.randn(*y_train_shape)

    model = TProcessRegressor(normalize_y=normalize_y)

    # FIXME: before fitting, the estimator does not have information regarding
    # the number of targets and default to 1. This is inconsistent with the shape
    # provided after `fit`. This assert should be made once the following issue
    # is fixed:
    # https://github.com/scikit-learn/scikit-learn/issues/22430
    # y_samples = model.sample_y(X_test, n_samples=n_samples_y_test)
    # assert y_samples.shape == y_test_shape

    model.fit(X_train, y_train)

    y_samples = model.sample_y(X_test, n_samples=n_samples_y_test)
    assert y_samples.shape == y_test_shape


@pytest.mark.parametrize("n_targets", [None, 1, 2, 3])
@pytest.mark.parametrize("n_samples", [1, 5])
def test_sample_y_shape_with_prior(n_targets, n_samples):
    """Check the output shape of `sample_y` is consistent before and after `fit`."""
    rng = np.random.RandomState(1024)

    X = rng.randn(10, 3)
    y = rng.randn(10, n_targets if n_targets is not None else 1)

    model = TProcessRegressor(n_targets=n_targets)
    shape_before_fit = model.sample_y(X, n_samples=n_samples).shape
    model.fit(X, y)
    shape_after_fit = model.sample_y(X, n_samples=n_samples).shape
    assert shape_before_fit == shape_after_fit


@pytest.mark.parametrize("n_targets", [None, 1, 2, 3])
def test_predict_shape_with_prior(n_targets):
    """Check the output shape of `predict` with prior distribution."""
    rng = np.random.RandomState(1024)

    n_sample = 10
    X = rng.randn(n_sample, 3)
    y = rng.randn(n_sample, n_targets if n_targets is not None else 1)

    model = TProcessRegressor(n_targets=n_targets)
    mean_prior, cov_prior = model.predict(X, return_cov=True)
    _, std_prior = model.predict(X, return_std=True)

    model.fit(X, y)
    mean_post, cov_post = model.predict(X, return_cov=True)
    _, std_post = model.predict(X, return_std=True)

    assert mean_prior.shape == mean_post.shape
    assert cov_prior.shape == cov_post.shape
    assert std_prior.shape == std_post.shape


def test_n_targets_error():
    """Check that an error is raised when the number of targets seen at fit is
    inconsistent with n_targets.
    """
    rng = np.random.RandomState(0)
    X = rng.randn(10, 3)
    y = rng.randn(10, 2)

    model = TProcessRegressor(n_targets=1)
    with pytest.raises(ValueError, match="The number of targets seen in `y`"):
        model.fit(X, y)


class CustomKernel(C):
    """
    A custom kernel that has a diag method that returns the first column of the
    input matrix X. This is a helper for the test to check that the input
    matrix X is not mutated.
    """

    def diag(self, X):
        return X[:, 0]


def test_tpr_predict_input_not_modified():
    """
    Check that the input X is not modified by the predict method of the
    TProcessRegressor when setting return_std=True.

    Non-regression test for:
    https://github.com/scikit-learn/scikit-learn/issues/24340
    """
    tpr = TProcessRegressor(kernel=CustomKernel()).fit(X, y)

    X2_copy = np.copy(X2)
    _, _ = tpr.predict(X2, return_std=True)

    assert_allclose(X2, X2_copy)


def test_degrees_of_freedom():
    """Checks that degrees of freedom are being added with each observation"""
    # Tests starting degrees of freedom
    kernel = RBF(length_scale=1.0)
    tpr_default = TProcessRegressor(kernel, optimizer=None)
    tpr_5 = TProcessRegressor(kernel, v=5, optimizer=None)
    assert tpr_default.v == 3
    assert tpr_5.v == 5

    # Tests degrees of freedom are being added to given observations
    tpr_default.fit(X, y)
    tpr_5.fit(X, y)
    assert tpr_default.v_n == 9
    assert tpr_5.v_n == 11


def test_multi_output_degrees_of_freedom():
    """
    Checks the degrees of freedom are beinf added correctly when using
    mulitple outputs
    """
    y_2d = np.vstack((y, y * 2)).T
    kernel = RBF(length_scale=1.0)

    tpr_1 = TProcessRegressor(kernel=kernel, optimizer=None, normalize_y=False)
    assert tpr_1.v == 3
    tpr_1.fit(X, y)
    assert tpr_1.v_n == 9


def test_measn_match_gp():
    """
    Checks the degrees of freedom are beinf added correctly when using
    mulitple outputs
    """

    kernel = RBF(length_scale=1.0, length_scale_bounds="fixed")
    tpr = TProcessRegressor(kernel=kernel).fit(X, y)
    gpr = TProcessRegressor(kernel=kernel).fit(X, y)

    mean_tp, _ = tpr.predict(X, return_std=True)
    mean_gp, _ = gpr.predict(X, return_std=True)

    assert_almost_equal(mean_tp, mean_gp)
