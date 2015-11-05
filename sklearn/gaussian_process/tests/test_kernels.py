"""Testing for kernels for Gaussian processes."""

# Author: Jan Hendrik Metzen <jhm@informatik.uni-bremen.de>
# Licence: BSD 3 clause

from collections import Hashable
from sklearn.externals.funcsigs import signature

import numpy as np

from sklearn.gaussian_process.kernels import _approx_fprime

from sklearn.metrics.pairwise \
    import PAIRWISE_KERNEL_FUNCTIONS, euclidean_distances, pairwise_kernels
from sklearn.gaussian_process.kernels \
    import (RBF, Matern, RationalQuadratic, ExpSineSquared, DotProduct,
            ConstantKernel, WhiteKernel, PairwiseKernel, KernelOperator,
            Exponentiation)
from sklearn.base import clone

from sklearn.utils.testing import (assert_equal, assert_almost_equal,
                                   assert_not_equal, assert_array_equal,
                                   assert_array_almost_equal)


X = np.random.RandomState(0).normal(0, 1, (5, 2))
Y = np.random.RandomState(0).normal(0, 1, (6, 2))

kernel_white = RBF(length_scale=2.0) + WhiteKernel(noise_level=3.0)
kernels = [RBF(length_scale=2.0), RBF(length_scale_bounds=(0.5, 2.0)),
           ConstantKernel(constant_value=10.0),
           2.0 * RBF(length_scale=0.33, length_scale_bounds="fixed"),
           2.0 * RBF(length_scale=0.5), kernel_white,
           2.0 * RBF(length_scale=[0.5, 2.0]),
           2.0 * Matern(length_scale=0.33, length_scale_bounds="fixed"),
           2.0 * Matern(length_scale=0.5, nu=0.5),
           2.0 * Matern(length_scale=1.5, nu=1.5),
           2.0 * Matern(length_scale=2.5, nu=2.5),
           2.0 * Matern(length_scale=[0.5, 2.0], nu=0.5),
           3.0 * Matern(length_scale=[2.0, 0.5], nu=1.5),
           4.0 * Matern(length_scale=[0.5, 0.5], nu=2.5),
           RationalQuadratic(length_scale=0.5, alpha=1.5),
           ExpSineSquared(length_scale=0.5, periodicity=1.5),
           DotProduct(sigma_0=2.0), DotProduct(sigma_0=2.0) ** 2]
for metric in PAIRWISE_KERNEL_FUNCTIONS:
    if metric in ["additive_chi2", "chi2"]:
        continue
    kernels.append(PairwiseKernel(gamma=1.0, metric=metric))


def test_kernel_gradient():
    """ Compare analytic and numeric gradient of kernels. """
    for kernel in kernels:
        K, K_gradient = kernel(X, eval_gradient=True)

        assert_equal(K_gradient.shape[0], X.shape[0])
        assert_equal(K_gradient.shape[1], X.shape[0])
        assert_equal(K_gradient.shape[2], kernel.theta.shape[0])

        def eval_kernel_for_theta(theta):
            kernel_clone = kernel.clone_with_theta(theta)
            K = kernel_clone(X, eval_gradient=False)
            return K

        K_gradient_approx = \
            _approx_fprime(kernel.theta, eval_kernel_for_theta, 1e-10)

        assert_almost_equal(K_gradient, K_gradient_approx, 4)


def test_kernel_theta():
    """ Check that parameter vector theta of kernel is set correctly. """
    for kernel in kernels:
        if isinstance(kernel, KernelOperator) \
           or isinstance(kernel, Exponentiation):  # skip non-basic kernels
            continue
        theta = kernel.theta
        _, K_gradient = kernel(X, eval_gradient=True)

        # Determine kernel parameters that contribute to theta
        init_sign = signature(kernel.__class__.__init__).parameters.values()
        args = [p.name for p in init_sign if p.name != 'self']
        theta_vars = map(lambda s: s.rstrip("_bounds"),
                         filter(lambda s: s.endswith("_bounds"), args))
        assert_equal(
            set(hyperparameter.name
                for hyperparameter in kernel.hyperparameters),
            set(theta_vars))

        # Check that values returned in theta are consistent with
        # hyperparameter values (being their logarithms)
        for i, hyperparameter in enumerate(kernel.hyperparameters):
            assert_equal(theta[i],
                         np.log(getattr(kernel, hyperparameter.name)))

        # Fixed kernel parameters must be excluded from theta and gradient.
        for i, hyperparameter in enumerate(kernel.hyperparameters):
            # create copy with certain hyperparameter fixed
            params = kernel.get_params()
            params[hyperparameter.name + "_bounds"] = "fixed"
            kernel_class = kernel.__class__
            new_kernel = kernel_class(**params)
            # Check that theta and K_gradient are identical with the fixed
            # dimension left out
            _, K_gradient_new = new_kernel(X, eval_gradient=True)
            assert_equal(theta.shape[0], new_kernel.theta.shape[0] + 1)
            assert_equal(K_gradient.shape[2], K_gradient_new.shape[2] + 1)
            if i > 0:
                assert_equal(theta[:i], new_kernel.theta[:i])
                assert_array_equal(K_gradient[..., :i],
                                   K_gradient_new[..., :i])
            if i + 1 < len(kernel.hyperparameters):
                assert_equal(theta[i+1:], new_kernel.theta[i:])
                assert_array_equal(K_gradient[..., i+1:],
                                   K_gradient_new[..., i:])

        # Check that values of theta are modified correctly
        for i, hyperparameter in enumerate(kernel.hyperparameters):
            theta[i] = np.log(42)
            kernel.theta = theta
            assert_almost_equal(getattr(kernel, hyperparameter.name), 42)

            setattr(kernel, hyperparameter.name, 43)
            assert_almost_equal(kernel.theta[i], np.log(43))


def test_auto_vs_cross():
    """ Auto-correlation and cross-correlation should be consistent. """
    for kernel in kernels:
        if kernel == kernel_white:
            continue  # Identity is not satisfied on diagonal
        K_auto = kernel(X)
        K_cross = kernel(X, X)
        assert_almost_equal(K_auto, K_cross, 5)


def test_kernel_diag():
    """ Test that diag method of kernel returns consistent results. """
    for kernel in kernels:
        K_call_diag = np.diag(kernel(X))
        K_diag = kernel.diag(X)
        assert_almost_equal(K_call_diag, K_diag, 5)


def test_kernel_operator_commutative():
    """ Adding kernels and multiplying kernels should be commutative. """
    # Check addition
    assert_almost_equal((RBF(2.0) + 1.0)(X),
                        (1.0 + RBF(2.0))(X))

    # Check multiplication
    assert_almost_equal((3.0 * RBF(2.0))(X),
                        (RBF(2.0) * 3.0)(X))


def test_kernel_anisotropic():
    """ Anisotropic kernel should be consistent with isotropic kernels."""
    kernel = 3.0 * RBF([0.5, 2.0])

    K = kernel(X)
    X1 = np.array(X)
    X1[:, 0] *= 4
    K1 = 3.0 * RBF(2.0)(X1)
    assert_almost_equal(K, K1)

    X2 = np.array(X)
    X2[:, 1] /= 4
    K2 = 3.0 * RBF(0.5)(X2)
    assert_almost_equal(K, K2)

    # Check getting and setting via theta
    kernel.theta = kernel.theta + np.log(2)
    assert_array_equal(kernel.theta, np.log([6.0, 1.0, 4.0]))
    assert_array_equal(kernel.k2.length_scale, [1.0, 4.0])


def test_kernel_stationary():
    """ Test stationarity of kernels."""
    for kernel in kernels:
        if not kernel.is_stationary():
            continue
        K = kernel(X, X + 1)
        assert_almost_equal(K[0, 0], np.diag(K))


def test_kernel_clone():
    """ Test that sklearn's clone works correctly on kernels. """
    for kernel in kernels:
        kernel_cloned = clone(kernel)

        assert_equal(kernel, kernel_cloned)
        assert_not_equal(id(kernel), id(kernel_cloned))
        for attr in kernel.__dict__.keys():
            attr_value = getattr(kernel, attr)
            attr_value_cloned = getattr(kernel_cloned, attr)
            if attr.startswith("hyperparameter_"):
                assert_equal(attr_value.name, attr_value_cloned.name)
                assert_equal(attr_value.value_type,
                             attr_value_cloned.value_type)
                assert_array_equal(attr_value.bounds,
                                   attr_value_cloned.bounds)
                assert_equal(attr_value.n_elements,
                             attr_value_cloned.n_elements)
            elif np.iterable(attr_value):
                for i in range(len(attr_value)):
                    if np.iterable(attr_value[i]):
                        assert_array_equal(attr_value[i],
                                           attr_value_cloned[i])
                    else:
                        assert_equal(attr_value[i], attr_value_cloned[i])
            else:
                assert_equal(attr_value, attr_value_cloned)
            if not isinstance(attr_value, Hashable):
                # modifiable attributes must not be identical
                assert_not_equal(id(attr_value), id(attr_value_cloned))


def test_matern_kernel():
    """ Test consistency of Matern kernel for special values of nu. """
    K = Matern(nu=1.5, length_scale=1.0)(X)
    # the diagonal elements of a matern kernel are 1
    assert_array_almost_equal(np.diag(K), np.ones(X.shape[0]))
    # matern kernel for coef0==0.5 is equal to absolute exponential kernel
    K_absexp = np.exp(-euclidean_distances(X, X, squared=False))
    K = Matern(nu=0.5, length_scale=1.0)(X)
    assert_array_almost_equal(K, K_absexp)
    # test that special cases of matern kernel (coef0 in [0.5, 1.5, 2.5])
    # result in nearly identical results as the general case for coef0 in
    # [0.5 + tiny, 1.5 + tiny, 2.5 + tiny]
    tiny = 1e-10
    for nu in [0.5, 1.5, 2.5]:
        K1 = Matern(nu=nu, length_scale=1.0)(X)
        K2 = Matern(nu=nu + tiny, length_scale=1.0)(X)
        assert_array_almost_equal(K1, K2)


def test_kernel_versus_pairwise():
    """Check that GP kernels can also be used as pairwise kernels."""
    for kernel in kernels:
        # Test auto-kernel
        if kernel != kernel_white:
            # For WhiteKernel: k(X) != k(X,X). This is assumed by
            # pairwise_kernels
            K1 = kernel(X)
            K2 = pairwise_kernels(X, metric=kernel)
            assert_array_almost_equal(K1, K2)

        # Test cross-kernel
        K1 = kernel(X, Y)
        K2 = pairwise_kernels(X, Y, metric=kernel)
        assert_array_almost_equal(K1, K2)


def test_set_get_params():
    """Check that set_params()/get_params() is consistent with kernel.theta."""
    for kernel in kernels:
        # Test get_params()
        index = 0
        params = kernel.get_params()
        for hyperparameter in kernel.hyperparameters:
            if hyperparameter.bounds is "fixed":
                continue
            size = hyperparameter.n_elements
            if size > 1:  # anisotropic kernels
                assert_almost_equal(np.exp(kernel.theta[index:index+size]),
                                    params[hyperparameter.name])
                index += size
            else:
                assert_almost_equal(np.exp(kernel.theta[index]),
                                    params[hyperparameter.name])
                index += 1
        # Test set_params()
        index = 0
        value = 10  # arbitrary value
        for hyperparameter in kernel.hyperparameters:
            if hyperparameter.bounds is "fixed":
                continue
            size = hyperparameter.n_elements
            if size > 1:  # anisotropic kernels
                kernel.set_params(**{hyperparameter.name: [value]*size})
                assert_almost_equal(np.exp(kernel.theta[index:index+size]),
                                    [value]*size)
                index += size
            else:
                kernel.set_params(**{hyperparameter.name: value})
                assert_almost_equal(np.exp(kernel.theta[index]), value)
                index += 1
