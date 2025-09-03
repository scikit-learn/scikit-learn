import numpy
from sklearn.tree import DecisionTreeRegressor
from sklearn.tree._utils import _py_precompute_absolute_errors

if False:
    def test_first_split():
        reg = DecisionTreeRegressor(max_depth=1, criterion='absolute_error')

        def mae_min(y, w):
            return min((np.abs(y - yi) * w).sum() for yi in y)

        def leaves_mae(l, y):
            return np.array([mae_min(y[l == i], w[l == i]) for i in np.unique(l)])

        X = np.array([
            [ 2.38, 3.13],
            [-0.87, 0.24],
            [ 3.42, 2.74],
            [ 1.43, 2.57],
            [ 0.86, 0.26]
        ])
        y = np.array([0.784, 0.654, 1.125, 2.010, 0.614])
        w = np.array([0.622, 1.356, 1.206, 0.912, 1.424])
        
        leaves = reg.fit(X, y, sample_weight=w).apply(X)
        print(leaves)
        assert leaves_mae(leaves, y).sum() < 1.1


def sample_X_y_w(n):
    x_true = (numpy.random.rand(n) > 0.5).astype(float)
    X = numpy.array([
        numpy.random.randn(n) + 2*x_true,
        numpy.round(2*numpy.random.rand(n) + 2*x_true, 2)
    ]).T
    X_pred = numpy.array([
        numpy.random.randn(n) + 2*x_true,
        2*numpy.random.rand(n) + 2*x_true
    ]).T
    y = numpy.random.rand(n) + (numpy.random.rand(n) + 0.5) * x_true
    w = 0.5 + numpy.random.rand(n)
    return X, y, w, X_pred



def test_absolute_errors_precomputation_function():
    """
    Test the main bit of logic of the MAE(RegressionCriterion) class 
    (used by DecisionTreeRegressor())

    The implemation of the criterion "repose" on an efficient precomputation
    of left/right children absolute error for each split. This test verifies this
    part of the computation, in case of major refactor of the MAE class, it can be safely removed
    """

    def compute_abs_error(y: numpy.ndarray, w: numpy.ndarray):
        # 1) compute the weighted median
        # i.e. once ordered by y, search for i such that:
        # sum(w[:i]) <= 1/2  and sum(w[i+1:]) <= 1/2
        sorter = numpy.argsort(y)
        wc = numpy.cumsum(w[sorter])
        idx = numpy.searchsorted(wc, wc[-1] / 2)
        median = y[sorter[idx]]
        print(y, median)
        # 2) compute the AE
        return (numpy.abs(y - median) * w).sum()

    def compute_prefix_abs_errors_naive(y: numpy.ndarray, w: numpy.ndarray):
        y = y.ravel()
        return numpy.array([compute_abs_error(y[:i], w[:i]) for i in range(1, y.size + 1)])


    for n in [3, 5, 10, 20, 100, 300]:
        y = numpy.random.uniform(size=(n, 1))
        w = numpy.random.rand(n)
        indices = numpy.arange(n)
        abs_errors = _py_precompute_absolute_errors(y, w, indices)
        expected = compute_prefix_abs_errors_naive(y, w)
        assert numpy.allclose(abs_errors, expected)

        abs_errors = _py_precompute_absolute_errors(y, w, indices, suffix=True)
        expected = compute_prefix_abs_errors_naive(y[::-1], w[::-1])[::-1]
        assert numpy.allclose(abs_errors, expected)

        x = numpy.random.rand(n)
        indices = numpy.argsort(x)
        w[:] = 1
        y_sorted = y[indices]
        w_sorted = w[indices]

        abs_errors = _py_precompute_absolute_errors(y, w, indices)
        expected = compute_prefix_abs_errors_naive(y_sorted, w_sorted)
        assert numpy.allclose(abs_errors, expected)

        abs_errors = _py_precompute_absolute_errors(y, w, indices, suffix=True)
        expected = compute_prefix_abs_errors_naive(y_sorted[::-1], w_sorted[::-1])[::-1]
        assert numpy.allclose(abs_errors, expected)



def test_first_split():

    def mae_min(y, w):
        return min((numpy.abs(y - yi) * w).sum() for yi in y)

    def leaves_mae(l, y, w=None):
        if w is None:
            w = numpy.ones(y.size)
        return numpy.array([mae_min(y[l == i], w[l == i]) for i in numpy.unique(l)])

    it = 0
    for n in [5]*100 + [10]*100 + [100]*100 + [1000]*10 + [10_000]*3:
        it += 1
        X, y, w, _ = sample_X_y_w(n)

        reg = DecisionTreeRegressor(max_depth=1, criterion='absolute_error')
        sk_leaves = reg.fit(X, y, sample_weight=w).apply(X)
        h_leaves = fit_apply(X, y, X_apply=X, sample_weights=w)
        are_leaves_the_same = (sk_leaves == h_leaves).all() or (sk_leaves == (3 - h_leaves)).all()
        if not are_leaves_the_same:
            sk_mae = leaves_mae(sk_leaves, y, w).sum()
            h_mae = leaves_mae(h_leaves, y, w).sum()
            assert numpy.isclose(sk_mae, h_mae), it

    for n in [5]*100 + [10]*100 + [100]*100 + [1000]*10 + [10_000]*3:
        it += 1
        X, y, _, _ = sample_X_y_w(n)
        reg = DecisionTreeRegressor(max_depth=1, criterion='absolute_error')
        sk_leaves = reg.fit(X, y).apply(X)
        h_leaves = fit_apply(X, y, X)
        are_leaves_the_same = (sk_leaves == h_leaves).all() or (sk_leaves == (3 - h_leaves)).all()
        if not are_leaves_the_same:
            sk_mae = leaves_mae(sk_leaves, y).sum()
            h_mae = leaves_mae(h_leaves, y).sum()
            assert numpy.isclose(sk_mae, h_mae), it


def fit_apply(X, y, X_apply=None, sample_weights=None):
    if sample_weights is None:
        sample_weights = numpy.ones(y.size)
    X_apply = X if X_apply is None else X_apply
    best_mae, best_feature, best_threshold = numpy.inf, -1, numpy.nan
    for k, x in enumerate(X.T):
        threshold, split_mae = min_mae_split(x, y, sample_weights)
        if split_mae < best_mae:
            best_mae = split_mae
            best_feature = k
            best_threshold = threshold
    return (X_apply[:, best_feature] >= best_threshold) + 1


def min_mae_split(x, y, w):
    """
    Find the best split of x that minimizes the sum of left and right MAEs.

    Sorts and deduplicates x, y, w, then computes the MAE for all possible splits using splits_left_mae.
    Returns the split value and the corresponding MAE.

    Parameters
    ----------
    x : np.ndarray
        Feature values.
    y : np.ndarray
        Target values.
    w : np.ndarray
        Sample weights.

    Returns
    -------
    x_split : float
        The value of x at which to split.
    split_mae : float
        The minimum sum of left and right MAEs.
    """
    sorter = numpy.argsort(x)
    x = x[sorter]
    y = y[sorter]
    w = w[sorter]
    prefix_maes = compute_prefix_maes(y, w)
    suffix_maes = compute_prefix_maes(y[::-1], w[::-1])[::-1]
    maes = prefix_maes + suffix_maes  # size: n-1
    maes[x[:-1] == x[1:]] = numpy.inf  # impossible to split between 2 points that are exactly equals
    best_split = numpy.argmin(maes)
    # Choose split point between best_split and its neighbor with lower MAE
    x_split = (x[best_split] + x[best_split + 1]) / 2
    split_mae = maes[best_split]
    return x_split, split_mae


def compute_prefix_maes(y: numpy.ndarray, w: numpy.ndarray):
    """
    Compute the minimum mean absolute error (MAE) for all (y[:i], w[:i]) with i ranging in [1, n-1]
    O(n log n) complexity, expect for patological cases (w growing faster than x^2)

    Parameters
    ----------
    y : numpy.ndarray
        Array of target values
    w : numpy.ndarray
        Array of sample weights.
    Returns
    -------
    maes : numpy.ndarray
        Prefix array of MAE values
    """
    n = y.size
    above = WeightedHeap(n, True)  # Min-heap for values above the median
    below = WeightedHeap(n, False)  # Max-heap for values below the median
    maes = numpy.full(n-1, numpy.inf)
    for i in range(n - 1):
        # Insert y[i] into the appropriate heap
        if above.empty() or y[i] > below.top():
            above.push(y[i], w[i])
        else:
            below.push(y[i], w[i])

        half_weight = (above.total_weight + below.total_weight) / 2
        # Rebalance the heaps, we want to ensure that:
        # above.total_weight >= 1/2 and above.total_weight - above.top_weight() <= 1/2
        # which ensures that above.top() is a weighted median of the heap
        # and in particular, an argmin for the MAE
        while above.total_weight < half_weight:
            yt, wt = below.pop()
            above.push(yt, wt)
        while above.total_weight - above.top_weight() > half_weight:
            yt, wt = above.pop()
            below.push(yt, wt)

        median = above.top()  # Current weighted median
        # Compute MAE for this split
        maes[i] = (
            (below.total_weight - above.total_weight) * median
            - below.weighted_sum
            + above.weighted_sum
        )
    return maes


class WeightedHeap:

    def __init__(self, max_size, min_heap=True):
        self.heap = numpy.zeros(max_size, dtype=numpy.float64)
        self.weights = numpy.zeros(max_size, dtype=numpy.float64)
        self.total_weight = 0
        self.weighted_sum = 0
        self.size = 0
        self.min_heap = min_heap

    def empty(self):
        return self.size == 0

    def push(self, val, weight):
        self.heap[self.size] = val if self.min_heap else -val
        self.weights[self.size] = weight
        self.total_weight += weight
        self.weighted_sum += val * weight
        self.size += 1
        self._perc_up(self.size - 1)

    def swap(self, i, j):
        self.heap[i], self.heap[j] = self.heap[j], self.heap[i]
        self.weights[i], self.weights[j] = self.weights[j], self.weights[i]

    def top(self):
        return self.heap[0] if self.min_heap else -self.heap[0]

    def top_weight(self):
        return self.weights[0]

    def pop(self):
        retv = self.top()
        retw = self.top_weight()
        self.size -= 1
        self.total_weight -= retw
        self.weighted_sum -= retv * retw
        self.heap[0] = self.heap[self.size]
        self.weights[0] = self.weights[self.size]
        self._perc_down(0)
        return retv, retw

    def _perc_up(self, i):
        p = (i - 1) >> 1
        while p >= 0:
            if self.heap[i] < self.heap[p]:
                self.swap(i, p)
            i = p
            p = (i - 1) >> 1

    def _perc_down(self, i):
        while (i << 1) + 2 <= self.size:
            mc_i = self._min_child_node(i)
            if self.heap[i] > self.heap[mc_i]:
                self.swap(i, mc_i)
            i = mc_i

    def _min_child_node(self, i):
        if (i << 1) + 2 == self.size:
            return (i << 1) | 1
        else:
            if self.heap[(i << 1) | 1] < self.heap[(i << 1) + 2]:
                return (i << 1) | 1
            else:
                return (i << 1) + 2

