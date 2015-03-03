# Author: Doug Coleman <doug.coleman@gmail.com>
# License: BSD 3 clause

from collections import defaultdict
import warnings
import numbers

import numpy as np
from ..utils.validation import check_random_state
from ..utils.random import sample_without_replacement

__all__ = ['resample_labels']


def _collect_indices(y):
    """Collects a list of indices for each element.
    Returns labels and indices where those labels appear in y.
    >>> _collect_indices(np.array([1, 2, 2, 3, 3, 3, 4, 4, 4, 4]))
    ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
    ((1, 2, 3, 4), ([0], [1, 2], [3, 4, 5], [6, 7, 8, 9]))
    """
    d = defaultdict(list)
    for k, v in enumerate(y):
        d[v].append(k)
    labels, indices = zip(*sorted(d.items()))
    return labels, indices


def _circular_sample(initial_n_samples, target_n_samples, random_state=None):
    """Sample without replacement from array in a loop.
    Take each sample once per loop except for the last loop, which takes
    (n % len(array)) samples at random. Results are returned unshuffled.
    >>> _circular_sample(3, 8, random_state=46)
    array([0, 0, 1, 1, 2, 2, 1, 0])
    """
    random_state = check_random_state(random_state)
    n_full_sets = int(target_n_samples // initial_n_samples)
    n_remainder = target_n_samples % initial_n_samples
    sample_indices = np.empty(target_n_samples, dtype=int)
    last_loop_index = n_full_sets * initial_n_samples
    sample_indices[:last_loop_index] =\
        np.repeat(np.arange(initial_n_samples), n_full_sets)
    if n_remainder > 0:
        sample_indices[last_loop_index:] = \
            sample_without_replacement(initial_n_samples, n_remainder,
                                       random_state=random_state)
    return sample_indices


def _fair_array_counts(n_samples, n_classes, random_state=None):
    """Tries to fairly partition n_samples between n_classes.
    If this cannot be done fairly, +1 is added `remainder` times
    to the counts for random arrays until a total of `n_samples` is
    reached.

    >>> _fair_array_counts(5, 3, random_state=43)
    array([2, 1, 2])
    """
    if n_classes > n_samples:
        raise ValueError("The number of classes is greater"
                         " than the number of samples requested")
    sample_size = n_samples // n_classes
    sample_size_rem = n_samples % n_classes
    counts = np.repeat(sample_size, n_classes)
    if sample_size_rem > 0:
        counts[:sample_size_rem] += 1
        # Shuffle so the class inbalance varies between runs
        random_state = check_random_state(random_state)
        random_state.shuffle(counts)
    return counts


def _scale_n_samples(scaling, n):
    """Helper function to scale the number of samples."""
    if scaling is None:
        return n
    else:
        if isinstance(scaling, numbers.Number) and scaling < 0:
            raise ValueError("Scaling must be nonnegative: %s" % scaling)
        elif isinstance(scaling, float):
            return scaling * n
        elif isinstance(scaling, (numbers.Integral, np.integer)):
            return scaling
        else:
            raise ValueError("Invalid value for scaling, must be "
                             "float, int, or None: %s" % scaling)


def _weighted_sample(probas, n_samples, random_state=None):
    """Select indices from n_samples with a weighted probability.

    Parameters
    ---------
    probas : array-like
        Array of probabilities summing to 1.
    n_samples : integer
        The number of samples to draw from at random.
    random_state : int, or RandomState instance (optional)
        Control the sampling for reproducible behavior.
    """
    random_state = check_random_state(random_state)
    if abs(np.sum(probas) - 1.0) > .011:
        raise ValueError("Label distribution probabilites must sum to 1")
    cum_probas = np.cumsum(probas)
    cum_probas[-1] = 1  # ensure that the probabilities sum to 1
    space = np.linspace(0, 1, 10000)
    weighted_indices = np.searchsorted(cum_probas, space)
    rs = random_state.randint(0, high=len(space), size=n_samples)
    return weighted_indices[rs]


def resample_labels(y, method=None, scaling=None, replace=False,
                    shuffle=False, random_state=None):
    """Resamples a classes array `y` and returns an array of indices

    The default behavior is to output the same ``y``. The additional
    parameters control the desired class distribution of the indices and the
    number of samples in the output.

    Parameters
    ----------
    y : array-like, shape (n_samples,)
        Target classes. Pass in the entire classes array so that this
        function can work on the class distribution.
    method : "balance", "oversample", "undersample", dict (optional)
        None outputs samples with the same class distribution as `y`.
        "balance" rebalances the classes to be equally distributed,
            over/undersampling for `len(y)` samples by default.
        "oversample" grows all classes to the count of the largest class.
        "undersample" shrinks all classes to the count of the smallest class.
        dict with pairs of class, probability with values summing to 1
    scaling : integer, float (optional)
        Number of samples to return.
        None outputs the same number of samples.
        `integer` is an absolute number of samples.
        `float` is a scale factor.
    replace : boolean (False by default)
        Sample with replacement when True.
    shuffle : boolean (False by default)
        Shuffle the indices before returning them. This option can add
        significant overhead, so it is disabled by default.
    random_state : int, or RandomState instance (optional)
        Control the sampling for reproducible behavior.

    Returns
    -------
    indices : array-like, shape (n_samples,)
        Indices sampled from the dataset respecting a class distribution
        controlled by this function's parameters.

    Examples
    --------
    Sample without replacement to reduce the size of a dataset by half
    and keep the same class distribution. Note how to apply the indices to X.

    >>> from sklearn.preprocessing import resample_labels
    >>> import numpy as np
    >>> X = np.array([[100], [120], [130], [110], [130], [110]])
    >>> y = np.array([10, 12, 13, 11, 13, 11])
    >>> indices = resample_labels(y, scaling=.5, random_state=333)
    >>> indices, X[indices], y[indices]
    ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
    (array([4, 5, 3]), array([[130], [110], [110]]), array([13, 11, 11]))

    Sample with replacement the dataset to 1.5 times its size and balance
    the class counts.

    >>> y = np.array([30, 30, 30, 10, 20, 30])
    >>> indices = resample_labels(y, method="balance", scaling=1.5,
    ...                           replace=True, random_state=335)
    >>> indices, y[indices]
    ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
    (array([3, 3, 3, 4, 4, 4, 2, 0, 1]),
     array([10, 10, 10, 20, 20, 20, 30, 30, 30]))

    Oversample all classes to the max class count of three samples each.

    >>> y = np.array([1, 2, 2, 3, 3, 3])
    >>> indices = resample_labels(y, method="oversample", random_state=333)
    >>> indices, y[indices]
    ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
    (array([0, 0, 0, 1, 2, 1, 3, 4, 5]),
     array([1, 1, 1, 2, 2, 2, 3, 3, 3]))

    Undersample all classes to the min class count of one sample each and also
    scale the number of samples by two.

    >>> y = np.array([1, 2, 2, 3, 3, 3])
    >>> indices = resample_labels(y, method="undersample", scaling=2.0,
    ...                           random_state=333)
    >>> indices, y[indices]
    (array([0, 0, 1, 2, 3, 4]), array([1, 1, 2, 2, 3, 3]))

    Sample twelve times with a probability dict.

    >>> y = np.array([1, 2, 3])
    >>> indices = resample_labels(y, method={1: .1, 2: .1, 3: .8},
    ...     scaling=12, random_state=337, shuffle=True)
    >>> indices, y[indices]
    ... # doctest: +NORMALIZE_WHITESPACE, +ELLIPSIS
    (array([2, 1, 2, 0, 0, 2, 1, 2, 2, 2, 2, 2]),
     array([3, 2, 3, 1, 1, 3, 2, 3, 3, 3, 3, 3]))
    """
    random_state = check_random_state(random_state)

    if method is None:
        n_samples = int(_scale_n_samples(scaling, len(y)))
        if replace:
            # already shuffled after this call
            sample_indices = random_state.randint(0, len(y), n_samples)
        else:
            sample_indices = _circular_sample(len(y), n_samples, random_state)
            if shuffle:
                random_state.shuffle(sample_indices)
        return sample_indices

    labels, indices = _collect_indices(y)

    if method in ('balance', 'oversample', 'undersample'):
        if method == 'balance':
            n_samples = _scale_n_samples(scaling, len(y))
        else:
            if method == 'oversample':
                count = max(len(a) for a in indices)
            else:
                count = min(len(a) for a in indices)
            n_samples = _scale_n_samples(scaling, count * len(labels))
        counts = _fair_array_counts(n_samples, len(labels), random_state)

    elif isinstance(method, dict):
        n_samples = int(_scale_n_samples(scaling, len(y)))
        proba = dict((k, v) for k, v in list(method.items()) if v > 0)
        desired_classes = np.asarray(list(proba.keys()))
        desired_probs = np.asarray(list(proba.values()))
        diff = set(desired_classes) - set(labels)
        if len(diff) > 0:
            raise ValueError("Can't make desired distribution: "
                             "some classes in `proba` dict are not in `y`: %s"
                             % list(diff))
        bincounts = np.bincount(_weighted_sample(desired_probs,
                                                n_samples,
                                                random_state))
        indices = [indices[k] for k, v in enumerate(bincounts) if v > 0]
        counts = bincounts[bincounts > 0]

    else:
        raise ValueError("Invalid value for method: %s" % method)

    if replace:
        sample_indices = \
            [[array[i] for i in
                random_state.randint(0, len(array), int(count))]
                for array, count in zip(indices, counts)]
    else:
        sample_indices = \
            [[array[i] for i in
                _circular_sample(len(array), int(count), random_state)]
                for array, count in zip(indices, counts)]
    sample_indices = np.concatenate(sample_indices)

    if shuffle:
        random_state.shuffle(sample_indices)
    return sample_indices
