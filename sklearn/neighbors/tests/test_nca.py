import sys
import numpy as np
from numpy.testing import assert_array_equal
from sklearn.exceptions import ConvergenceWarning
from sklearn.externals.six import StringIO
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_raises, assert_equal, \
    assert_raise_message, assert_warns_message, assert_true
from sklearn.datasets import load_iris, make_classification
from sklearn.neighbors.nca import NeighborhoodComponentsAnalysis
from sklearn.metrics import pairwise_distances


rng = check_random_state(0)
# load and shuffle iris dataset
iris = load_iris()
perm = rng.permutation(iris.target.size)
iris_data = iris.data[perm]
iris_target = iris.target[perm]
EPS = np.finfo(float).eps


def test_simple_example():
    """Test on a simple example.

    Puts four points in the input space where the opposite labels points are
    next to each other. After transform the same labels points should be next
    to each other.

    """
    X = np.array([[0, 0], [0, 1], [2, 0], [2, 1]])
    y = np.array([1, 0, 1, 0])
    nca = NeighborhoodComponentsAnalysis(n_components=2, init='identity',
                                         random_state=42)
    nca.fit(X, y)
    Xansformed = nca.transform(X)
    np.testing.assert_equal(pairwise_distances(Xansformed).argsort()[:, 1],
                            np.array([2, 3, 0, 1]))


def test_finite_differences():
    r"""Test gradient of loss function

    Test if the gradient is correct by computing the relative difference
    between the projected gradient PG:

    .. math::

        PG = \mathbf d^{\top} \cdot \nabla
        \mathcal L(\mathbf x)

    and the finite differences FD:

    .. math::

        FD = \frac{\mathcal L(\mathbf x + \epsilon \mathbf d) -
        \mathcal L(\mathbf x - \epsilon \mathbf d)}{2 \epsilon}


    where :math:`d` is a random direction (random vector of shape `n_features`,
    and norm 1), :math:`\epsilon` is a very small number, :math:`\mathcal L` is
    the loss function and :math:`\nabla \mathcal L` is its gradient. This
    relative difference should be zero:

    .. math ::

        \frac{|PG -FD|}{|PG|} = 0


    """
    # Initialize `transformation`, `X` and `y` and `NCA`
    X = iris_data
    y = iris_target
    point = rng.randn(rng.randint(1, X.shape[1] + 1), X.shape[1])
    nca = NeighborhoodComponentsAnalysis(init=point)

    X, y, init = nca._validate_params(X, y)
    mask = y[:, np.newaxis] == y[np.newaxis, :]  # (n_samples, n_samples)
    nca.n_iter_ = 0

    point = nca._initialize(X, init)
    # compute the gradient at `point`
    _, gradient = nca._loss_grad_lbfgs(point, X, mask)

    # create a random direction of norm 1
    random_direction = rng.randn(*point.shape)
    random_direction /= np.linalg.norm(random_direction)

    # computes projected gradient
    projected_gradient = random_direction.ravel().dot(
                                      gradient.ravel())

    # compute finite differences
    eps = 1e-5
    right_loss, _ = nca._loss_grad_lbfgs(point + eps * random_direction, X,
                                         mask)
    left_loss, _ = nca._loss_grad_lbfgs(point - eps * random_direction, X,
                                        mask)
    finite_differences = 1 / (2 * eps) * (right_loss - left_loss)

    # compute relative error
    relative_error = np.abs(finite_differences - projected_gradient) / \
        np.abs(projected_gradient)
    np.testing.assert_almost_equal(relative_error, 0.)


def test_params_validation():
    # Test that invalid parameters raise value error
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]
    NCA = NeighborhoodComponentsAnalysis

    # TypeError
    assert_raises(TypeError, NCA(max_iter='21').fit, X, y)
    assert_raises(TypeError, NCA(verbose='true').fit, X, y)
    assert_raises(TypeError, NCA(tol=1).fit, X, y)
    assert_raises(TypeError, NCA(n_components='invalid').fit, X, y)
    assert_raises(TypeError, NCA(warm_start=1).fit, X, y)

    # ValueError
    assert_raise_message(ValueError,
                         "`init` must be 'pca', 'identity', 'random' or a "
                         "numpy array of shape (n_components, n_features).",
                         NCA(init=1).fit, X, y)
    assert_raise_message(ValueError,
                         '`max_iter`= -1, must be >= 1.',
                         NCA(max_iter=-1).fit, X, y)

    init = np.random.rand(5, 3)
    assert_raise_message(ValueError,
                         'The output dimensionality ({}) of the given linear '
                         'transformation `init` cannot be greater than its '
                         'input dimensionality ({}).'
                         .format(init.shape[0], init.shape[1]),
                         NCA(init=init).fit, X, y)

    n_components = 10
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) cannot be greater '
                         'than the given data dimensionality ({})!'
                         .format(n_components, X.shape[1]),
                         NCA(n_components=n_components).fit, X, y)


def test_transformation_dimensions():
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]

    # Fail if transformation input dimension does not match inputs dimensions
    transformation = np.array([[1, 2], [3, 4]])
    assert_raises(ValueError,
                  NeighborhoodComponentsAnalysis(init=transformation).fit,
                  X, y)

    # Fail if transformation output dimension is larger than
    # transformation input dimension
    transformation = np.array([[1, 2], [3, 4], [5, 6]])
    # len(transformation) > len(transformation[0])
    assert_raises(ValueError,
                  NeighborhoodComponentsAnalysis(init=transformation).fit,
                  X, y)

    # Pass otherwise
    transformation = np.arange(9).reshape(3, 3)
    NeighborhoodComponentsAnalysis(init=transformation).fit(X, y)


def test_n_components():
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]

    init = np.random.rand(X.shape[1] - 1, 3)

    # n_components = X.shape[1] != transformation.shape[0]
    n_components = X.shape[1]
    nca = NeighborhoodComponentsAnalysis(init=init, n_components=n_components)
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) does not match '
                         'the output dimensionality of the given '
                         'linear transformation `init` ({})!'
                         .format(n_components, init.shape[0]),
                         nca.fit, X, y)

    # n_components > X.shape[1]
    n_components = X.shape[1] + 2
    nca = NeighborhoodComponentsAnalysis(init=init, n_components=n_components)
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) cannot be greater '
                         'than the given data dimensionality ({})!'
                         .format(n_components, X.shape[1]),
                         nca.fit, X, y)

    # n_components < X.shape[1]
    nca = NeighborhoodComponentsAnalysis(n_components=2, init='identity')
    nca.fit(X, y)


def test_init_transformation():
    X, y = make_classification(n_samples=30, n_features=5,
                               n_redundant=0, random_state=0)

    # Start learning from scratch
    nca = NeighborhoodComponentsAnalysis(init='identity')
    nca.fit(X, y)

    # Initialize with random
    nca_random = NeighborhoodComponentsAnalysis(init='random')
    nca_random.fit(X, y)

    # Initialize with PCA
    nca_pca = NeighborhoodComponentsAnalysis(init='pca')
    nca_pca.fit(X, y)

    init = np.random.rand(X.shape[1], X.shape[1])
    nca = NeighborhoodComponentsAnalysis(init=init)
    nca.fit(X, y)

    # init.shape[1] must match X.shape[1]
    init = np.random.rand(X.shape[1], X.shape[1] + 1)
    nca = NeighborhoodComponentsAnalysis(init=init)
    assert_raise_message(ValueError,
                         'The input dimensionality ({}) of the given '
                         'linear transformation `init` must match the '
                         'dimensionality of the given inputs `X` ({}).'
                         .format(init.shape[1], X.shape[1]),
                         nca.fit, X, y)

    # init.shape[0] must be <= init.shape[1]
    init = np.random.rand(X.shape[1] + 1, X.shape[1])
    nca = NeighborhoodComponentsAnalysis(init=init)
    assert_raise_message(ValueError,
                         'The output dimensionality ({}) of the given '
                         'linear transformation `init` cannot be '
                         'greater than its input dimensionality ({}).'
                         .format(init.shape[0], init.shape[1]),
                         nca.fit, X, y)

    # init.shape[0] must match n_components
    init = np.random.rand(X.shape[1], X.shape[1])
    n_components = X.shape[1] - 2
    nca = NeighborhoodComponentsAnalysis(init=init, n_components=n_components)
    assert_raise_message(ValueError,
                         'The preferred embedding dimensionality '
                         '`n_components` ({}) does not match '
                         'the output dimensionality of the given '
                         'linear transformation `init` ({})!'
                         .format(n_components, init.shape[0]),
                         nca.fit, X, y)


def test_warm_start_validation():
    X, y = make_classification(n_samples=30, n_features=5, n_classes=4,
                               n_redundant=0, n_informative=5, random_state=0)

    nca = NeighborhoodComponentsAnalysis(warm_start=True, max_iter=5)
    nca.fit(X, y)

    X_less_features, y = \
        make_classification(n_samples=30, n_features=4, n_classes=4,
                            n_redundant=0, n_informative=4, random_state=0)
    assert_raise_message(ValueError,
                         'The new inputs dimensionality ({}) does not '
                         'match the input dimensionality of the '
                         'previously learned transformation ({}).'
                         .format(X_less_features.shape[1],
                                 nca.components_.shape[1]),
                         nca.fit, X_less_features, y)


def test_warm_start_effectiveness():
    # A 1-iteration second fit on same data should give almost same result
    # with warm starting, and quite different result without warm starting.

    X, y = make_classification(n_samples=30, n_features=5,
                               n_redundant=0, random_state=0)
    n_iter = 10

    nca_warm = NeighborhoodComponentsAnalysis(warm_start=True,
                                              max_iter=n_iter, random_state=0)
    nca_warm.fit(X, y)
    transformation_warm = nca_warm.components_
    nca_warm.max_iter = 1
    nca_warm.fit(X, y)
    transformation_warm_plus_one = nca_warm.components_

    nca_cold = NeighborhoodComponentsAnalysis(warm_start=False,
                                              max_iter=n_iter, random_state=0)
    nca_cold.fit(X, y)
    transformation_cold = nca_cold.components_
    nca_cold.max_iter = 1
    nca_cold.fit(X, y)
    transformation_cold_plus_one = nca_cold.components_

    diff_warm = np.sum(np.abs(transformation_warm_plus_one -
                              transformation_warm))
    diff_cold = np.sum(np.abs(transformation_cold_plus_one -
                              transformation_cold))

    assert_true(diff_warm < 2.0,
                "Transformer changed significantly after one iteration even "
                "though it was warm-started.")

    assert_true(diff_cold > diff_warm,
                "Cold-started transformer changed less significantly than "
                "warm-started transformer after one iteration.")


def test_verbose():
    # assert there is proper output when verbose = 1
    old_stdout = sys.stdout
    sys.stdout = StringIO()

    nca = NeighborhoodComponentsAnalysis(verbose=1)
    try:
        nca.fit(iris_data, iris_target)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    # check output
    assert("[NeighborhoodComponentsAnalysis]" in out)
    assert("Finding principal components" in out)
    assert ("Finding principal components" in out)
    assert ("Training took" in out)

    # assert by default there is no output (verbose=0)
    old_stdout = sys.stdout
    sys.stdout = StringIO()

    nca = NeighborhoodComponentsAnalysis()
    try:
        nca.fit(iris_data, iris_target)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    # check output
    assert(out == '')


def test_singleton_class():
    X = iris_data
    y = iris_target

    # one singleton class
    singleton_class = 1
    ind_singleton, = np.where(y == singleton_class)
    y[ind_singleton] = 2
    y[ind_singleton[0]] = singleton_class

    nca = NeighborhoodComponentsAnalysis(max_iter=30)
    nca.fit(X, y)

    # One non-singleton class
    ind_1, = np.where(y == 1)
    ind_2, = np.where(y == 2)
    y[ind_1] = 0
    y[ind_1[0]] = 1
    y[ind_2] = 0
    y[ind_2[0]] = 2

    nca = NeighborhoodComponentsAnalysis(max_iter=30)
    nca.fit(X, y)

    # Only singleton classes
    ind_0, = np.where(y == 0)
    ind_1, = np.where(y == 1)
    ind_2, = np.where(y == 2)
    X = X[[ind_0[0], ind_1[0], ind_2[0]]]
    y = y[[ind_0[0], ind_1[0], ind_2[0]]]

    nca = NeighborhoodComponentsAnalysis(init='identity', max_iter=30)
    nca.fit(X, y)
    assert_array_equal(X, nca.transform(X))


def test_one_class():
    X = iris_data[iris_target == 0]
    y = iris_target[iris_target == 0]

    nca = NeighborhoodComponentsAnalysis(max_iter=30,
                                         n_components=X.shape[1],
                                         init='identity')
    nca.fit(X, y)
    assert_array_equal(X, nca.transform(X))


def test_callback():
    X = iris_data
    y = iris_target

    nca = NeighborhoodComponentsAnalysis(callback='my_cb')
    assert_raises(ValueError, nca.fit, X, y)

    max_iter = 10

    def my_cb(transformation, n_iter):
        rem_iter = max_iter - n_iter
        print('{} iterations remaining...'.format(rem_iter))

    # assert that my_cb is called
    old_stdout = sys.stdout
    sys.stdout = StringIO()

    nca = NeighborhoodComponentsAnalysis(max_iter=max_iter,
                                         callback=my_cb, verbose=1)
    try:
        nca.fit(iris_data, iris_target)
    finally:
        out = sys.stdout.getvalue()
        sys.stdout.close()
        sys.stdout = old_stdout

    # check output
    assert('{} iterations remaining...'.format(max_iter-1) in out)


def test_store_opt_result():
    X = iris_data
    y = iris_target

    nca = NeighborhoodComponentsAnalysis(max_iter=5,
                                         store_opt_result=True)
    nca.fit(X, y)
    transformation = nca.opt_result_.x
    assert_equal(transformation.size, X.shape[1]**2)


def test_convergence_warning():

    nca = NeighborhoodComponentsAnalysis(max_iter=2, verbose=1)
    cls_name = nca.__class__.__name__
    assert_warns_message(ConvergenceWarning,
                         '[{}] NCA did not converge'.format(cls_name),
                         nca.fit, iris_data, iris_target)
