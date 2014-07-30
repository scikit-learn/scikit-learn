import numpy as np

from sklearn.isotonic import check_increasing, isotonic_regression,\
    IsotonicRegression

from sklearn.utils.testing import assert_raises, assert_array_equal,\
    assert_true, assert_false, assert_equal

from sklearn.utils.testing import assert_warns_message, assert_no_warnings


def test_check_increasing_up():
    x = [0, 1, 2, 3, 4, 5]
    y = [0, 1.5, 2.77, 8.99, 8.99, 50]

    # Check that we got increasing=True and no warnings
    is_increasing = assert_no_warnings(check_increasing, x, y)
    assert_true(is_increasing)


def test_check_increasing_up_extreme():
    x = [0, 1, 2, 3, 4, 5]
    y = [0, 1, 2, 3, 4, 5]

    # Check that we got increasing=True and no warnings
    is_increasing = assert_no_warnings(check_increasing, x, y)
    assert_true(is_increasing)


def test_check_increasing_down():
    x = [0, 1, 2, 3, 4, 5]
    y = [0, -1.5, -2.77, -8.99, -8.99, -50]

    # Check that we got increasing=False and no warnings
    is_increasing = assert_no_warnings(check_increasing, x, y)
    assert_false(is_increasing)


def test_check_increasing_down_extreme():
    x = [0, 1, 2, 3, 4, 5]
    y = [0, -1, -2, -3, -4, -5]

    # Check that we got increasing=False and no warnings
    is_increasing = assert_no_warnings(check_increasing, x, y)
    assert_false(is_increasing)


def test_check_ci_warn():
    x = [0, 1, 2, 3, 4, 5]
    y = [0, -1, 2, -3, 4, -5]

    # Check that we got increasing=False and CI interval warning
    is_increasing = assert_warns_message(UserWarning, "interval",
                                         check_increasing,
                                         x, y)

    assert_false(is_increasing)


def test_isotonic_regression():
    y = np.array([3, 7, 5, 9, 8, 7, 10])
    y_ = np.array([3, 6, 6, 8, 8, 8, 10])
    assert_array_equal(y_, isotonic_regression(y))

    x = np.arange(len(y))
    ir = IsotonicRegression(y_min=0., y_max=1.)
    ir.fit(x, y)
    assert_array_equal(ir.fit(x, y).transform(x), ir.fit_transform(x, y))
    assert_array_equal(ir.transform(x), ir.predict(x))

    # check that it is immune to permutation
    perm = np.random.permutation(len(y))
    ir = IsotonicRegression(y_min=0., y_max=1.)
    assert_array_equal(ir.fit_transform(x[perm], y[perm]),
                       ir.fit_transform(x, y)[perm])
    assert_array_equal(ir.transform(x[perm]), ir.transform(x)[perm])

    # check it doesn't change y when all x are equal:
    ir = IsotonicRegression()
    assert_array_equal(ir.fit_transform(np.ones(len(x)), y), y)


def test_isotonic_regression_reversed():
    y = np.array([10, 9, 10, 7, 6, 6.1, 5])
    y_ = IsotonicRegression(increasing=False).fit_transform(
        np.arange(len(y)), y)
    assert_array_equal(np.ones(y_[:-1].shape), ((y_[:-1] - y_[1:]) >= 0))


def test_isotonic_regression_auto_decreasing():
    # Set y and x for decreasing
    y = np.array([10, 9, 10, 7, 6, 6.1, 5])
    x = np.arange(len(y))

    # Create model and fit_transform
    ir = IsotonicRegression(increasing='auto')
    y_ = assert_no_warnings(ir.fit_transform, x, y)

    # Check that relationship decreases
    is_increasing = y_[0] < y_[-1]
    assert_false(is_increasing)


def test_isotonic_regression_auto_increasing():
    # Set y and x for decreasing
    y = np.array([5, 6.1, 6, 7, 10, 9, 10])
    x = np.arange(len(y))

    # Create model and fit_transform
    ir = IsotonicRegression(increasing='auto')
    y_ = assert_no_warnings(ir.fit_transform, x, y)

    # Check that relationship increases
    is_increasing = y_[0] < y_[-1]
    assert_true(is_increasing)


def test_assert_raises_exceptions():
    ir = IsotonicRegression()
    rng = np.random.RandomState(42)
    assert_raises(ValueError, ir.fit, [0, 1, 2], [5, 7, 3], [0.1, 0.6])
    assert_raises(ValueError, ir.fit, [0, 1, 2], [5, 7])
    assert_raises(ValueError, ir.fit, rng.randn(3, 10), [0, 1, 2])
    assert_raises(ValueError, ir.transform, rng.randn(3, 10))


def test_isotonic_sample_weight_parameter_default_value():
    # check if default value of sample_weight parameter is one
    ir = IsotonicRegression()
    # random test data
    rng = np.random.RandomState(42)
    n = 100
    x = np.arange(n)
    y = rng.randint(-50, 50, size=(n,)) + 50. * np.log(1 + np.arange(n))
    # check if value is correctly used
    weights = np.ones(n)
    y_set_value = ir.fit_transform(x, y, sample_weight=weights)
    y_default_value = ir.fit_transform(x, y)

    assert_array_equal(y_set_value, y_default_value)


def test_isotonic_min_max_boundaries():
    # check if min value is used correctly
    ir = IsotonicRegression(y_min=2, y_max=4)
    n = 6
    x = np.arange(n)
    y = np.arange(n)
    y_test = [2, 2, 2, 3, 4, 4]
    y_result = np.round(ir.fit_transform(x, y))
    assert_array_equal(y_result, y_test)


def test_isotonic_sample_weight():
    ir = IsotonicRegression()
    x = [1, 2, 3, 4, 5, 6, 7]
    y = [1, 41, 51, 1, 2, 5, 24]
    sample_weight = [1, 2, 3, 4, 5, 6, 7]
    expected_y = [1, 13.95, 13.95, 13.95, 13.95, 13.95, 24]
    received_y = ir.fit_transform(x, y, sample_weight=sample_weight)

    assert_array_equal(expected_y, received_y)


def test_isotonic_regression_oob_raise():
    # Set y and x
    y = np.array([3, 7, 5, 9, 8, 7, 10])
    x = np.arange(len(y))

    # Create model and fit
    ir = IsotonicRegression(increasing='auto', out_of_bounds="raise")
    ir.fit(x, y)

    # Check that an exception is thrown
    assert_raises(ValueError, ir.predict, [min(x)-10, max(x)+10])


def test_isotonic_regression_oob_clip():
    # Set y and x
    y = np.array([3, 7, 5, 9, 8, 7, 10])
    x = np.arange(len(y))

    # Create model and fit
    ir = IsotonicRegression(increasing='auto', out_of_bounds="clip")
    ir.fit(x, y)

    # Predict from  training and test x and check that min/max match.
    y1 = ir.predict([min(x) - 10, max(x) + 10])
    y2 = ir.predict(x)
    assert_equal(max(y1), max(y2))
    assert_equal(min(y1), min(y2))


def test_isotonic_regression_oob_nan():
    # Set y and x
    y = np.array([3, 7, 5, 9, 8, 7, 10])
    x = np.arange(len(y))

    # Create model and fit
    ir = IsotonicRegression(increasing='auto', out_of_bounds="nan")
    ir.fit(x, y)

    # Predict from  training and test x and check that we have two NaNs.
    y1 = ir.predict([min(x)-10, max(x)+10])
    assert_equal(sum(np.isnan(y1)), 2)


def test_isotonic_regression_oob_bad():
    # Set y and x
    y = np.array([3, 7, 5, 9, 8, 7, 10])
    x = np.arange(len(y))

    # Create model and fit
    ir = IsotonicRegression(increasing='auto', out_of_bounds="xyz")

    # Make sure that we throw an error for bad out_of_bounds value
    assert_raises(ValueError, ir.fit, x, y)


def test_isotonic_regression_oob_bad_after():
    # Set y and x
    y = np.array([3, 7, 5, 9, 8, 7, 10])
    x = np.arange(len(y))

    # Create model and fit
    ir = IsotonicRegression(increasing='auto', out_of_bounds="raise")

    # Make sure that we throw an error for bad out_of_bounds value in transform
    ir.fit(x, y)
    ir.out_of_bounds = "xyz"
    assert_raises(ValueError, ir.transform, x)


if __name__ == "__main__":
    import nose
    nose.run(argv=['', __file__])
