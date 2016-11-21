import numpy as np

from sklearn.datasets import ValueDropper
from sklearn.datasets import make_classification, make_regression
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_true
from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal
from sklearn.utils.testing import assert_raise_message
from sklearn.utils.testing import assert_raises_regexp
from sklearn.preprocessing import LabelEncoder


def test_value_dropper_mnar_clf():
    # Test drop probabilites when missing distribution is
    # given for classification problems
    n_samples, n_features = 1000, 5
    n_values = n_samples * n_features
    X, y = make_classification(n_samples=n_samples,
                               n_classes=4,
                               n_features=n_features,
                               n_informative=5,
                               n_redundant=0,
                               n_repeated=0,
                               random_state=0)
    le = LabelEncoder().fit(['a', 'z', 'b', 'j'])
    y_str = le.inverse_transform(y)
    y_int = y

    for y in (y_int, y_str):
        classes = np.unique(y)

        # Inplace dropping of values

        # Samples from class 0 will have a drop probability of 0.1
        vd = ValueDropper(missing_distribution={classes[0]: 0.1},
                          missing_values=np.nan, random_state=0)
        X_dropped = vd.transform(X, y)

        # Check the total drop fraction
        assert_almost_equal(np.isnan(X_dropped).ravel().sum() /
                            float(n_values), 0.1)

        # All the missing values are from y == 0
        assert_almost_equal(np.isnan(
            X_dropped[y == classes[0]]).ravel().sum() / float(n_values), 0.1)

        # and no missing values from y != 0
        assert_almost_equal(np.isnan(
            X_dropped[y != classes[0]]).ravel().sum() / float(n_values), 0.)

        # Samples from class 0 will have a drop probabilty of 0.3
        # but spread unevenly across features as given by the
        # list of probabilities
        # And samples from class 1 will have a drop probability of 0.01
        # across all features

        missing_distribution = {classes[0]: [0.01, 0.005, 0.005, 0, 0],
                                classes[1]: 0.05}
        vd = ValueDropper(missing_distribution=missing_distribution,
                          missing_values=np.nan,
                          random_state=0)
        X_dropped = vd.transform(X, y)

        missing_mask = np.isnan(X_dropped)

        # Check that there are no missing values when y != {0 or 1}
        assert_equal(missing_mask[(y == classes[2])].ravel().sum(), 0)
        assert_equal(missing_mask[(y == classes[3])].ravel().sum(), 0)

        # Check that the drop probabilites when class == 1 is 0.1
        # across all features
        assert_equal(missing_mask[y == classes[1]].ravel().sum() /
                     float(n_values), 0.05)

        # Check that the drop probabilites when class == 0 is 2.1
        # across all features sum(missing_distribution[0])
        assert_equal(missing_mask[y == classes[0]].ravel().sum() /
                     float(n_values), 0.02)

        # Check that the features indexed 3 and 4 have no missing values
        # for class 0
        assert_equal(missing_mask[y == classes[0]][3, 4].ravel().sum(), 0.)

        # Check that feature indexed 0 has drop prob 0.2
        assert_equal(missing_mask[np.where(y == classes[0])[0],
                                  (0,)].ravel().sum() /
                     float(n_values), 0.01)

        # Check that feature indexed 1 and 2 both have drop prob 0.05
        # Check that feature indexed 0 has drop prob 0.2
        assert_equal(missing_mask[np.where(y == classes[0])[0],
                                  (1,)].ravel().sum() /
                     float(n_values), 0.005)
        assert_equal(missing_mask[np.where(y == classes[0])[0],
                                  (2,)].ravel().sum() /
                     float(n_values), 0.005)

        # Ensure scaling the missing_distribution by a factor of 2
        missing_distribution = {classes[0]: [0.02, 0.01, 0.01, 0, 0],
                                classes[1]: 0.1}
        vd = ValueDropper(missing_distribution=missing_distribution,
                          missing_values=np.nan, random_state=0)
        X_dropped2 = vd.transform(X, y)
        assert_true(np.all(np.isnan(X_dropped2[np.isnan(X_dropped)])))


def test_value_dropper_mnar_reg_error():
    X, y = make_regression(n_samples=10, random_state=0)

    assert_raise_message(ValueError,
                         "only for single target which is discrete"
                         " (classification tasks). The given target (y) is of "
                         "type continuous",
                         ValueDropper(missing_distribution={0: 0.2}).transform,
                         X, y)


def check_value_dropper_mcar(X, y):
    X_copy = X.copy()
    X_copy2 = X.copy()
    n_samples, n_features = X.shape
    n_values = n_samples * n_features

    # Inplace dropping of values; 0 correlation case.
    # For even indexed features missing probability is 0.03 and
    # for odd indexed ones 0.01
    missing_distribution = np.array([0.03, 0.01] * 5)
    vd = ValueDropper(missing_distribution=missing_distribution,
                      copy=False, random_state=0)
    vd.transform(X_copy, y)
    missing_mask = np.isnan(X_copy)

    global_missing_rate = missing_distribution.sum()

    # Check the global missing rate
    assert_almost_equal(missing_mask.ravel().sum() / float(n_values),
                        global_missing_rate)

    # Check the rate for all even indexed features
    assert_almost_equal(missing_mask[:, missing_distribution == 0.03]
                        .ravel().sum() / float(n_values),
                        0.03 * 5)

    # Check the rate for one even indexed feature
    assert_almost_equal(missing_mask[:, 0]
                        .ravel().sum() / float(n_values), 0.03)

    # Check the rate for all odd features
    assert_almost_equal(missing_mask[:, missing_distribution == 0.03]
                        .ravel().sum() / float(n_values),
                        0.03 * 5)

    # Check the rate for one odd indexed feature
    assert_almost_equal(missing_mask[:, 1]
                        .ravel().sum() / float(n_values), 0.01)

    # Let us drop 0.3 more fraction of values. This time not inplace
    # copy=True must be default
    vd = ValueDropper(missing_distribution=0.6, random_state=0)
    X_more_dropped = vd.transform(X_copy2, y)
    new_missing_mask = np.isnan(X_more_dropped)

    # Ensure X is not modified
    assert_array_almost_equal(X_copy2, X)

    # Ensure all the missing positions that were in the previous step also
    # exist when missing_distribution is scaled up
    # (Important for reproducibility)
    assert_true(np.all(new_missing_mask[missing_mask]))


def test_value_dropper_mcar():
    # Test missing fractions for MCAR case in a classification problem
    n_samples, n_features = 1000, 10
    X, y_int = make_classification(n_samples=n_samples,
                                   n_features=n_features, random_state=0)
    le = LabelEncoder().fit(['a', 'z'])
    y_str = le.inverse_transform(y_int)
    for y in (y_str, y_int):
        check_value_dropper_mcar(X, y)

    # Test missing fractions for MCAR case in a regression problem
    n_samples, n_features = 1000, 10
    X, y = make_regression(n_samples=n_samples, n_features=n_features,
                           random_state=0)
    check_value_dropper_mcar(X, y)


def test_value_dropper_errors():
    n_samples, n_features = 1000, 10
    X, y = make_classification(n_samples=n_samples,
                               n_classes=4,
                               n_features=n_features,
                               n_informative=5,
                               n_redundant=0,
                               n_repeated=0,
                               random_state=0)

    # Raise sensible error when sum of all probabilites in missing_distribution
    # exceeds or equals 1
    missing_distributions = (
        # NMAR cases
        {0: 0.25, 1: 0.25, 2: 0.25, 3: 0.25}, {0: 2., },
        {0: [0, 0, 0, 0, 0.24, 0, 0, 0, 0, 0.01], 1: 0.25, 2: [0.025, ] * 10,
         3: 0.25}, {0: [0.1] * 10, }, {0: 0.26, 1: [0.09, ] * 10},
        # MCAR cases
        [0, 0, 0, 0.2, 0.3, 0.1, 0, 0, 0, 0.5], 2.5, 1.5,
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0])
    for missing_distribution in missing_distributions:
        assert_raise_message(ValueError, "should sum up to less than 1",
                             ValueDropper(
                                 missing_distribution=missing_distribution)
                             .transform, X, y)

    # Accept only float values 0, 1, 2 all are incorrect values even as float
    # hence no point in accepting any int values
    assert_raise_message(ValueError,
                         "missing_distribution must be a float or  1D vector",
                         ValueDropper(missing_distribution=2).transform, X, y)
    assert_raise_message(ValueError,
                         "should either be a single float or an array",
                         ValueDropper(missing_distribution={0: 2})
                         .transform, X, y)

    wrong_missing_distributions_err_pairs = (
        # 1D vector with fewer or more than n_feature elements
        ([0.01, ] * 9, "does not conform to the number of features, 10"),
        ([0.01, ] * 11, "does not conform to the number of features, 10"),

        # Dict with labels having fewer or more than n_feature elements
        ({1: [0.01, ] * 9, },
         "for label, 1, does not conform to the number of features, 10"),

        ({0: [0.01, ] * 10, 1: [0.01, ] * 11},
         "for label, 1, does not conform to the number of features, 10"),

        # Dict having labels not present in y labels
        ({0: 0.025, 1: [0.0025, ] * 10, 2: 0.025, 3: 0.025, 4: 0.025},
         "y contains new labels: \[4\]"),

        # Incorrect dict or incorrect value
        ({0: {1: 0.2}, },
         "either be a single float or an array of shape \(n_features,\). "
         "\{1: 0.2.*\} was passed for class label 0"),

        ("foobar",
         "must be a float or  1D vector \(list, tuple or np.ndarray\)"
         " of shape \(n_features,\) or dict"))

    missing_distribution = {0: 0.025, 1: 0.025, 2: 0.025, 3: 0.025}
    for missing_distribution, err_msg in wrong_missing_distributions_err_pairs:
        assert_raises_regexp(ValueError, err_msg,
                             ValueDropper(
                                 missing_distribution=missing_distribution)
                             .transform, X, y)

    # When missing_distribution is a dict, but y is not given
    assert_raise_message(ValueError, "",
                         ValueDropper(
                             missing_distribution=missing_distribution)
                         .transform, X)
