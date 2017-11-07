import numpy as np
from numpy.testing import assert_array_equal
from sklearn.preprocessing import OneHotEncoder
from sklearn.utils import check_random_state
from sklearn.utils.testing import assert_raises, assert_equal
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
    nca = NeighborhoodComponentsAnalysis(n_features_out=2, init='identity',
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
    masks = OneHotEncoder(sparse=False,
                          dtype=bool).fit_transform(y[:, np.newaxis])
    nca.n_iter_ = 0

    point = nca._initialize(X, init)
    # compute the gradient at `point`
    _, gradient = nca._loss_grad_lbfgs(point, X, y, masks)

    # create a random direction of norm 1
    random_direction = rng.randn(*point.shape)
    random_direction /= np.linalg.norm(random_direction)

    # computes projected gradient
    projected_gradient = random_direction.ravel().dot(
                                      gradient.ravel())

    # compute finite differences
    eps = 1e-5
    right_loss, _ = nca._loss_grad_lbfgs(point + eps * random_direction, X, y,
                                         masks)
    left_loss, _ = nca._loss_grad_lbfgs(point - eps * random_direction, X, y,
                                        masks)
    finite_differences = 1/(2*eps) * (right_loss - left_loss)

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
    assert_raises(TypeError, NCA(n_features_out='invalid').fit, X, y)

    # ValueError
    assert_raises(ValueError, NCA(init=1).fit, X, y)
    assert_raises(ValueError, NCA(max_iter=-1).fit, X, y)

    fit_func = NCA(init=np.random.rand(5, 3)).fit
    assert_raises(ValueError, fit_func, X, y)
    assert_raises(ValueError, NCA(n_features_out=10).fit, X, y)


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


def test_n_features_out():
    X = np.arange(12).reshape(4, 3)
    y = [1, 1, 2, 2]

    transformation = np.array([[1, 2, 3], [4, 5, 6]])

    # n_features_out = X.shape[1] != transformation.shape[0]
    nca = NeighborhoodComponentsAnalysis(n_features_out=3, init=transformation)
    assert_raises(ValueError, nca.fit, X, y)

    # n_features_out > X.shape[1]
    nca = NeighborhoodComponentsAnalysis(n_features_out=5, init=transformation)
    assert_raises(ValueError, nca.fit, X, y)

    # n_features_out < X.shape[1]
    nca = NeighborhoodComponentsAnalysis(n_features_out=2, init='identity')
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
    assert_raises(ValueError, nca.fit, X, y)

    # init.shape[0] must be <= init.shape[1]
    init = np.random.rand(X.shape[1] + 1, X.shape[1])
    nca = NeighborhoodComponentsAnalysis(init=init)
    assert_raises(ValueError, nca.fit, X, y)

    # init.shape[0] must match n_features_out
    init = np.random.rand(X.shape[1], X.shape[1])
    nca = NeighborhoodComponentsAnalysis(n_features_out=X.shape[1] - 2,
                                         init=init)
    assert_raises(ValueError, nca.fit, X, y)


def test_verbose():
    nca = NeighborhoodComponentsAnalysis(verbose=1)
    nca.fit(iris_data, iris_target)
    # TODO: rather assert that some message is printed


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
                                         n_features_out=X.shape[1],
                                         init='identity')
    nca.fit(X, y)
    assert_array_equal(X, nca.transform(X))


def test_callable():
    X = iris_data
    y = iris_target

    nca = NeighborhoodComponentsAnalysis(callback='my_cb')
    assert_raises(ValueError, nca.fit, X, y)

    max_iter = 10

    def my_cb(transformation, n_iter):
        rem_iter = max_iter - n_iter
        print('{} iterations remaining...'.format(rem_iter))

    nca = NeighborhoodComponentsAnalysis(max_iter=max_iter,
                                         callback=my_cb, verbose=1)
    nca.fit(X, y)
    # TODO: rather assert that message is printed


def test_store_opt_result():
    X = iris_data
    y = iris_target

    nca = NeighborhoodComponentsAnalysis(max_iter=5,
                                         store_opt_result=True)
    nca.fit(X, y)
    transformation = nca.opt_result_.x
    assert_equal(transformation.size, X.shape[1]**2)
