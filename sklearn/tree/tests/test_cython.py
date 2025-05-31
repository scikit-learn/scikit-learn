import numpy as np

from sklearn.preprocessing import OneHotEncoder
from sklearn.tree import (
    DecisionTreeClassifier,
    ExtraTreeClassifier,
)
from sklearn.tree._partitioner import (
    DensePartitioner,
    py_breiman_sort_categories,
    py_init_node_split,
)


def test_breiman_sort_categories():
    # Example: 6 samples, 1 categorical feature with 4 possible categories (0, 1, 2, 3)
    X = np.array([[0], [1], [2], [1], [3], [2]], dtype=np.float32)
    y = np.array([[100], [20], [30], [40], [50], [60]], dtype=np.float64)
    sample_weight = np.array([1, 2, 1, 1, 1, 1], dtype=np.float64)
    n_categories = np.array([4], dtype=np.int32)
    missing_mask = np.zeros(1, dtype=np.uint8)
    samples = np.arange(X.shape[0], dtype=np.intp)
    feature_values = X[:, 0].copy()

    # Create the partitioner
    part = DensePartitioner(
        X, samples, feature_values, missing_mask, n_categories, breiman_shortcut=True
    )
    py_init_node_split(part, 0, X.shape[0])

    # Prepare cat_offset and sorted_cat arrays
    ncat_present = 4  # all categories present
    cat_offset = np.zeros(ncat_present, dtype=np.intp)
    sorted_cat = np.zeros(ncat_present, dtype=np.intp)

    # Call the method
    py_breiman_sort_categories(
        part, 0, 4, ncat_present, cat_offset, sorted_cat, y, sample_weight
    )

    # Compute expected means
    means = []
    for cat in range(4):
        mask = X[:, 0] == cat
        if np.sum(mask) == 0:
            means.append(np.inf)
        else:
            means.append(np.average(y[mask, 0], weights=sample_weight[mask]))
    expected_sorted = np.argsort(means)

    # Check that sorted_cat matches expected_sorted
    assert np.all(
        sorted_cat == expected_sorted
    ), f"Expected {expected_sorted}, got {sorted_cat}"
    print("Test passed! Sorted categories:", sorted_cat)


def test_breiman_sort_categories_with_cat_offset():
    # 10 possible categories, but only 2, 4, 7, 9 are present
    X = np.array([[2], [4], [7], [2], [9], [4], [7], [9]], dtype=np.float32)
    y = np.array([[10], [200], [30], [40], [50], [60], [70], [80]], dtype=np.float64)
    sample_weight = np.array([1, 2, 1, 1, 1, 1, 1, 1], dtype=np.float64)
    n_categories = np.array([10], dtype=np.int32)
    missing_mask = np.zeros(1, dtype=np.uint8)
    samples = np.arange(X.shape[0], dtype=np.intp)
    feature_values = X[:, 0].copy()

    # Create the partitioner
    part = DensePartitioner(
        X, samples, feature_values, missing_mask, n_categories, breiman_shortcut=True
    )
    py_init_node_split(part, 0, X.shape[0])

    # Only 4 categories are present: 2, 4, 7, 9
    present_cats = np.array([2, 4, 7, 9], dtype=np.intp)
    ncat_present = len(present_cats)
    # cat_offset maps local index to global category index
    cat_offset = present_cats.copy() - np.arange(ncat_present)
    sorted_cat = np.zeros(ncat_present, dtype=np.intp)

    # Call the method
    py_breiman_sort_categories(
        part, 0, 10, ncat_present, cat_offset, sorted_cat, y, sample_weight
    )

    # Compute expected means for present categories only
    means = []
    for cat in present_cats:
        mask = X[:, 0] == cat
        means.append(np.average(y[mask, 0], weights=sample_weight[mask]))
    expected_sorted_local = np.argsort(means)
    expected_sorted_global = present_cats[expected_sorted_local]

    # Check that sorted_cat matches expected_sorted_global
    assert np.all(
        sorted_cat == expected_sorted_global
    ), f"Expected {expected_sorted_global}, got {sorted_cat}"
    print("Test passed! Sorted categories (global indices):", sorted_cat)


def test_categorical_split_vs_onehot_tree_depth():
    # Simulate data: 8 categories, 1000 samples, binary classification,
    # non-ordinal mapping
    rng = np.random.RandomState(42)
    n_samples = 1000
    n_categories = 8
    X = rng.randint(0, n_categories, size=(n_samples, 1))
    # Shuffle categories and assign half to class 0, half to class 1
    categories = np.arange(n_categories)
    rng.shuffle(categories)
    class0_cats = set(categories[:4])
    class1_cats = set(categories[4:])
    y = np.array([0 if x[0] in class0_cats else 1 for x in X])

    # Train with categorical feature
    tree_cat = DecisionTreeClassifier(random_state=0, categorical_features=[0])
    tree_cat.fit(X, y)

    # Train with one-hot encoding
    ohe = OneHotEncoder(sparse_output=False, categories="auto")
    X_ohe = ohe.fit_transform(X)
    tree_ohe = DecisionTreeClassifier(random_state=0)
    tree_ohe.fit(X_ohe, y)

    # The categorical split should yield a shallower tree
    assert (
        tree_cat.get_depth() == 0
    ), f"Categorical split depth should be 1, got {tree_cat.get_depth()}"
    assert (
        tree_ohe.get_depth() >= 3
    ), f"One-hot tree depth should be at least 3, got {tree_ohe.get_depth()}"


CLF_TREES = {
    "DecisionTreeClassifier": DecisionTreeClassifier,
    "ExtraTreeClassifier": ExtraTreeClassifier,
}
# toy sample
X = [[-2, -1], [-1, -1], [-1, -2], [1, 1], [1, 2], [2, 1]]
y = [-1, -1, -1, 1, 1, 1]
T = [[-1, -1], [2, 2], [3, 2]]
true_result = [-1, 1, 1]
from sklearn.utils._testing import (
    assert_array_equal,
)


def test_weighted_classification_toy():
    # Check classification on a weighted toy dataset.
    for name, Tree in CLF_TREES.items():
        clf = Tree(random_state=0)

        clf.fit(X, y, sample_weight=np.ones(len(X)))
        assert_array_equal(clf.predict(T), true_result, "Failed with {0}".format(name))

        clf.fit(X, y, sample_weight=np.full(len(X), 0.5))
        assert_array_equal(clf.predict(T), true_result, "Failed with {0}".format(name))


test_weighted_classification_toy()
test_categorical_split_vs_onehot_tree_depth()
