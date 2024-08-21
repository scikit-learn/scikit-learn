from sklearn.datasets import load_iris
from sklearn.feature_selection._relieff import ReliefF
from sklearn.utils._testing import assert_allclose


def test_relieff_with_iris():
    """Test that ReliefF produces the expected outputs using the Iris dataset.
    The values are compared against those that are presented in the examples
    within the Matlab implementation.
    """
    iris_data = load_iris()
    X = iris_data.data
    y = iris_data.target

    expected_weights = [0.1399, 0.1226, 0.3590, 0.3754]
    relf = ReliefF(k=10, neighbors_algorithm="kd_tree")
    relf.fit(X, y)
    assert_allclose(relf.weights_, expected_weights, rtol=1e-3)
