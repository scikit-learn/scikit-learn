import numpy as np

from sklearn.utils.testing import assert_almost_equal
from sklearn.utils.testing import assert_array_equal
from sklearn.utils.testing import assert_equal
from sklearn.utils.testing import assert_raises
from sklearn.preprocessing.resample import _collect_indices
from sklearn.preprocessing.resample import _fair_array_counts
from sklearn.preprocessing.resample import resample_labels

y = np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4])


def test_collect_indices():
    labels, indices = _collect_indices(y)

    expected_labels = [1, 2, 3, 4]
    expected_indices = [[0], [1, 2], [3, 4, 5], [6, 7, 8, 9]]
    assert_array_equal(labels, expected_labels)
    assert_array_equal(indices, expected_indices)


def test_resample_labels_same_distribution():
    indices = resample_labels(y)
    assert_equal(len(indices), 10)
    assert_array_equal(y, y[np.sort(indices)])

    indices = resample_labels(y, replace=False)
    assert_equal(len(indices), 10)
    assert_array_equal(y, y[np.sort(indices)])

    indices = resample_labels(y, replace=True)
    assert_equal(len(indices), 10)

    indices = resample_labels(y, scaling=2.0, replace=False)
    assert_array_equal(y[np.sort(indices)], np.repeat(y, 2))

    indices = resample_labels(y, scaling=20, replace=False)
    assert_array_equal(y[np.sort(indices)], np.repeat(y, 2))

    indices = resample_labels(y, scaling=2.0, replace=False, shuffle=True)
    assert_array_equal(y[np.sort(indices)], np.sort(np.repeat(y, 2)))

    indices = resample_labels(y, scaling=20, replace=False, shuffle=True)
    assert_array_equal(y[np.sort(indices)], np.sort(np.repeat(y, 2)))

    indices = resample_labels(y, scaling=2.0, replace=True)
    assert_equal(len(indices), 20)

    indices = resample_labels(y, scaling=20, replace=True)
    assert_equal(len(indices), 20)


def test_resample_balanced():
    indices = resample_labels(y, scaling=2.0, method="balance",
                              replace=False)
    assert_equal(len(indices), 2 * len(y))

    indices = resample_labels(y, scaling=20, method="balance",
                              replace=False)
    assert_equal(len(indices), 2 * len(y))

    indices = resample_labels(y, scaling=2.0, method="balance",
                              replace=True)
    assert_equal(len(indices), 2 * len(y))

    indices = resample_labels(y, scaling=20, method="balance",
                              replace=True)
    assert_equal(len(indices), 2 * len(y))

    # [0,0,0,0, 1,1,1,1, 2,2,2,2, ...]
    z = np.concatenate([[i] * (4) for i in range(101)])
    check_mean_std(z, 50.0, 29.1547)

    # Should be equal because we sample everything and don't replace
    indices = resample_labels(z, method=None, replace=False,
                              random_state=42)
    assert_array_equal(z, z[indices])

    # Sample with replacement, should get nearly same mean/std as y
    assert_almost_equal(np.mean(z), 50.0, 3)
    indices = resample_labels(z, method=None, replace=True,
                              random_state=42)
    assert_almost_equal(np.mean(z[indices]), 51.2425, 3)
    assert_almost_equal(np.std(z[indices]), 28.5119, 3)

    indices = resample_labels(z, method=None, scaling=20000,
                              replace=True, random_state=42)
    assert_almost_equal(np.mean(z[indices]), 50.0067, 3)
    assert_almost_equal(np.std(z[indices]), 29.2240, 3)


def test_resample_labels_oversample():
    indices = resample_labels(y, method="oversample", replace=False)
    assert_equal(len(indices), 16)

    indices = resample_labels(y, method="oversample", replace=True)
    assert_equal(len(indices), 16)


def test_resample_labels_undersample():
    indices = resample_labels(y, method="undersample", replace=False)
    assert_equal(len(indices), 4)

    indices = resample_labels(y, method="undersample", replace=True)
    assert_equal(len(indices), 4)


def test_resample_labels_dict():
    indices = resample_labels(y, scaling=2.0,
                              method={1: .3, 2: .1, 3: .5, 4: .1},
                              random_state=43)
    assert_equal(len(indices), 2 * len(y))
    assert_array_equal(np.bincount(y[indices]), [0, 5, 1, 9, 5])

    indices = resample_labels(y, scaling=100,
                              method={1: .3, 2: .1, 3: .5, 4: .1},
                              random_state=42)
    assert_equal(len(indices), 100)
    assert_array_equal(np.bincount(y[indices]), [0, 28,  9, 54,  9])

    indices = resample_labels(y, scaling=100,
                              method={1: .3, 2: .7, 999: 0},
                              random_state=42)
    assert_equal(len(indices), 100)
    assert_array_equal(np.bincount(y[indices]), [0, 28, 72])

    indices = resample_labels(y, scaling=2.0,
                              method={1: .3, 2: .1, 3: .5, 4: .1},
                              replace=True, random_state=43)
    assert_equal(len(indices), 2 * len(y))
    assert_array_equal(np.bincount(y[indices]), [0, 5, 1, 9, 5])

    indices = resample_labels(y, scaling=100,
                              method={1: .3, 2: .1, 3: .5, 4: .1},
                              replace=True, random_state=42)
    assert_equal(len(indices), 100)
    assert_array_equal(np.bincount(y[indices]), [0, 28, 9, 54, 9])

    indices = resample_labels(y, scaling=100,
                              method={1: .3, 2: .7, 999: 0},
                              replace=True, random_state=42)
    assert_equal(len(indices), 100)
    assert_array_equal(np.bincount(y[indices]), [0, 28, 72])

    indices = resample_labels(y, scaling=100000,
                              method={1: .3, 2: .7, 999: 0},
                              replace=True,
                              random_state=42)
    assert_equal(len(indices), 100000)
    assert_array_equal(np.bincount(y[indices]), [0, 29972, 70028])


def test_resample_labels_invalid_parameters():
    assert_raises(ValueError, resample_labels, y, scaling=2.0,
                  method="badstring")
    assert_raises(ValueError, resample_labels, y, scaling=-2.0)
    assert_raises(ValueError, resample_labels, y, scaling=-2)
    assert_raises(ValueError, resample_labels, y, method="badstring")
    assert_raises(ValueError, resample_labels, y, method={1: .5})
    assert_raises(ValueError, resample_labels, y, method={})
    assert_raises(ValueError, resample_labels, y, method=555)
    assert_raises(ValueError, _fair_array_counts, 2, 5)
    assert_raises(ValueError, resample_labels, y, scaling="failme")
    assert_raises(ValueError, resample_labels,
                  y, method={1: .1, 2: .1, 30000: .8},
                  scaling=12, random_state=337, shuffle=True)


def check_mean_std(y, expected_mean, expected_std):
    indices = resample_labels(y, method="balance", replace=False,
                              random_state=42)
    assert_almost_equal(np.mean(y[indices]), expected_mean, 3)
    assert_almost_equal(np.std(y[indices]), expected_std, 3)

    indices = resample_labels(y, method="balance", replace=True,
                              random_state=42)
    assert_almost_equal(np.mean(y[indices]), expected_mean, 3)
    assert_almost_equal(np.std(y[indices]), expected_std, 3)

    # Scale the dataset and check the invariant
    indices = resample_labels(y, method="balance",
                              scaling=4.0, replace=False, random_state=42)
    assert_almost_equal(np.mean(y[indices]), expected_mean, 3)
    assert_almost_equal(np.std(y[indices]), expected_std, 3)

    indices = resample_labels(y, method="balance", scaling=4.0,
                              replace=True, random_state=42)
    assert_almost_equal(np.mean(y[indices]), expected_mean, 3)
    assert_almost_equal(np.std(y[indices]), expected_std, 3)


def test_resample_labels_linearly_unbalanced():
    # [0, 1,1, 2,2,2, ...]
    y = np.concatenate([[i] * (i + 1) for i in range(101)])
    check_mean_std(y, 50.0, 29.1547)

    # Should be equal because we sample everything and don't replace
    indices = resample_labels(y, method=None, replace=False,
                              random_state=42)
    assert_array_equal(y, y[indices])

    # Sample with replacement, should get nearly same mean/std as y
    assert_almost_equal(np.mean(y), 66.6666, 3)
    indices = resample_labels(y, method=None, replace=True,
                              random_state=42)
    assert_almost_equal(np.mean(y[indices]), 67.1574, 3)
    assert_almost_equal(np.std(y[indices]), 23.5837, 3)

    indices = resample_labels(y, method=None, scaling=10000,
                              replace=True, random_state=42)
    assert_almost_equal(np.mean(y[indices]), 67.0377, 3)
    assert_almost_equal(np.std(y[indices]), 23.6514, 3)


def test_resample_labels_quadratically_unbalanced():
    # [0, 1,1, 2,2,2,2, ...]
    y = np.concatenate([[i] * 2 ** i for i in range(11)])
    # Multiple of 11 for exact mean
    y = y[:2002]
    check_mean_std(y, 5.0, 3.1622)

    # Should be equal because we sample everything and don't replace
    indices = resample_labels(y, method=None, replace=False,
                              random_state=42)
    assert_array_equal(y, y[indices])

    # Sample with replacement, should get nearly same mean/std as y
    assert_almost_equal(np.mean(y), 8.9830, 3)
    indices = resample_labels(y, method=None, replace=True,
                              random_state=42)
    assert_almost_equal(np.mean(y[indices]), 8.9755, 3)
    assert_almost_equal(np.std(y[indices]), 1.4740, 3)

    indices = resample_labels(y, scaling=4004, method=None, replace=True,
                              random_state=42)
    assert_almost_equal(np.mean(y[indices]), 8.9822, 3)
    assert_almost_equal(np.std(y[indices]), 1.4315, 3)


def test_resample_labels_shuffled():
    assert_array_equal(
        resample_labels(y, replace=False, shuffle=False, random_state=42),
        np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9]))

    assert_array_equal(
        resample_labels(y, replace=False, shuffle=True, random_state=42),
        np.array([8, 1, 5, 0, 7, 2, 9, 4, 3, 6]))

    assert_array_equal(
        resample_labels(y, replace=True, shuffle=False, random_state=42),
        np.array([6, 3, 7, 4, 6, 9, 2, 6, 7, 4]))

    assert_array_equal(
        resample_labels(y, replace=True, shuffle=True, random_state=42),
        np.array([6, 3, 7, 4, 6, 9, 2, 6, 7, 4]))


def test_resample_labels_distribution_shuffled():
    assert_array_equal(
        resample_labels(y, method="balance", replace=False,
                        shuffle=False,
                        random_state=42),
        np.array([0, 0, 0, 1, 2, 3, 4, 5, 7, 9])
    )

    assert_array_equal(
        resample_labels(y, method="balance", replace=True,
                        shuffle=False,
                        random_state=42),
        np.array([0, 0, 0, 1, 2, 3, 3, 5, 7, 8])
    )

    assert_array_equal(
        resample_labels(y, method="balance", replace=False,
                        shuffle=True,
                        random_state=42),
        np.array([0, 0, 3, 1, 2, 5, 7, 9, 0, 4])
    )

    assert_array_equal(
        resample_labels(y, method="balance", replace=True,
                        shuffle=True,
                        random_state=42),
        np.array([0, 3, 8, 0, 7, 0, 1, 2, 5, 3])
    )


def test_resample_labels_dict_shuffled():
    assert_array_equal(
        resample_labels(y, method={1: .5, 2: .5},
                        replace=False, shuffle=False,
                        random_state=42),
        np.array([0, 0, 0, 1, 1, 1, 2, 2, 2, 1])
    )

    assert_array_equal(
        resample_labels(y, method={1: .5, 2: .5},
                        replace=False, shuffle=True,
                        random_state=42),
        np.array([1, 0, 0, 2, 2, 1, 1, 2, 0, 1])
    )

    assert_array_equal(
        resample_labels(y, method={1: .5, 2: .5},
                        replace=True, shuffle=False,
                        random_state=42),
        np.array([0, 0, 0, 2, 1, 2, 2, 2, 2, 2])
    )

    assert_array_equal(
        resample_labels(y, method={1: .5, 2: .5},
                        replace=True, shuffle=True,
                        random_state=42),
        np.array([2, 0, 2, 2, 2, 2, 2, 0, 1, 0])
    )
