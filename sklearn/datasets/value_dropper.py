# Author : Raghav RV <rvraghav93@gmail.com>
#
# Licence : BSD 3 clause

import numpy as np
import numbers

from sklearn.utils import check_array
from sklearn.utils import check_random_state
from sklearn.utils.multiclass import type_of_target

from sklearn.base import TransformerMixin
from sklearn.preprocessing import LabelEncoder


__all__ = ["ValueDropper"]


class ValueDropper(TransformerMixin):
    """Artificially insert NMAR or MCAR missing values into data.

    Where,

    NMAR/MNAR - Not Missing At Random / Missing Not At Random
        When the missingness is correlated with the class classes in the
        target (y) (and hence informative).

    MCAR - Missing Completely At Random
        When the missingness is completely random (and hence uninformative).

    If the missing type is NMAR, a ``missing_distribution``
    parameter can be passed to drop values conforming to that distribution.


    Parameters
    ----------

    missing_values : {"NaN" (or np.nan) | int}, default "NaN"
        The value to insert to indicate missingness.


    missing_distribution : dict of floats or dict of vector of floats
        If ``missing_distribution`` is a float within range [0, 1),
        it represents the absolute fraction of values that will be missing::
            0.2

        There will be ``0.2 * n_samples * n_features`` numbers of missing
        values in the data.

        Alternatively this refers to the probability of a value being dropped
        after transform. The values are dropped (approximately) uniformly
        across all labels and features. This type of missingness is referred
        to as MCAR.

        To vary the distribution across features or to prevent a feature from
        having missing values, individual probabilities for each feature can
        be specified as a 1D vector of shape ``(n_features,)``::
            [0.2, 0.2, 0]

        For the above example, the probability that a sample will have missing
        value in feature 0 is 0.2. In other words, after calling
        ``transform``, there are ``0.2 * n_samples`` randomly selected samples
        with value missing in feature 0 and 1 but no samples with
        missing values in feature 2.

        If missingness is not MCAR, ``missing_distribution`` can be used
        to specify the multinomial distribution of the newly dropped values
        across labels (and if needed across features) as given below.

        If ``missing_distribution`` is a dict of floats::
            {1: 0.02, 2: 0.03, 3: 0.05}

        The probability that a sample from class 1 will have a missing value
        is 0.02.

        The missing values are evenly spread across all the features.

        In other words, there are ``int(0.02 / n_features * n_samples)``
        randomly chosen samples of class 1 having missing values in each
        feature.

        If there are fewer than ``int(0.02 / n_features * n_samples)`` numbers
        of samples in class 1, an error is raised.

        Hence the total missing rate is ``0.02 + 0.03 + 0.05 = 0.1``.

        If ``missing_distribution`` is a dict of vectors (and scalars)::

            {0: 0.1,
             3: [0.1, 0.15, 0.15]}

        Note that the shape of the vector must be ``(n_features,)``

        There are ``0.1 / n_features * n_samples`` randomly chosen samples,
        for each feature, of class 1 having a missing value.

        There are 0 samples of class 1 and 2 having missing value in any
        feature.

        And There are ``0.1 * n_samples`` randomly chosen samples having
        missing value in feature 0, ``0.15 * n_samples`` randomly chosen
        samples having missing value in feature 1 and feature 2, each.

        A few notes:

        Note that the samples are randomly chosen "*for each feature*".

        The global ``missing_rate`` (fraction of missing values in the entire
        dataset) is calculated as the sum of all the individual probabilities.

        At the end of transform, the dataset will contain a total of
        ``(X.shape[0] * X.shape[1]) * missing_rate`` numbers of missing values.

    copy : bool, default False
        Whether to copy the data or work inplace.

    random_state : int, optional
        The seed for the numpy's random number generator.

        If ``random_state`` is set to an integer, the ``missing_distribution``
        can be scaled uniformly with the assumption that all the values
        dropped with a smaller scale will exist in the larger scaled version::
            missing_distribution_25pc = {0: 0.05, 3: [0.05, 0.075, 0.075]}
            missing_distribution_50pc = {0: 0.1, 3: [0.1, 0.15, 0.15]}

        The missing values dropped with ``missing_distribution_25pc`` will also
        exist in ``missing_distribution_50pc``.

        This guarantee does not apply when relative probabilities
        within the ``missing_distribution`` change between two settings::
            missing_distribution_25pc = {0: 0.05, 3: [0.05, 0.075, 0.075]}
            # The below dist. is not a scaled version of above
            missing_distribution_50pc = {0: 0.09, 3: [0.11, 0.15, 0.15]}


    Examples
    --------

    >>> import numpy as np
    >>> X = np.array([[0., 1., 2.],
    ...               [3., 4., 5.],
    ...               [6., 7., 8.],
    ...               [9., 0., 1.],
    ...               [2., 3., 4.],
    ...               [8., 9., 8.],
    ...               [1., 0., 5.],
    ...               [7., 8., 9.],
    ...               [5., 4., 3.],
    ...               [2., 1., 1.],
    ...               [1., 2., 3.]])
    >>> y = [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]
    >>> abs_missing_rate = 0.1
    >>> # Drop values across features 0 and 1 in
    >>> # the ratio of 4:1 from samples with class label as 1
    >>> # The final fraction values that will be missing
    >>> missing_distribution = {1: [0.8 * abs_missing_rate,
    ...                             0.2 * abs_missing_rate,
    ...                             0]}
    >>> vd = ValueDropper(missing_distribution=missing_distribution,
    ...                   random_state=0)
    >>> vd.transform(X, y)
    array([[  0.,   1.,   2.],
           [  3.,   4.,   5.],
           [  6.,   7.,   8.],
           [  9.,   0.,   1.],
           [  2.,   3.,   4.],
           [  8.,   9.,   8.],
           [  1.,   0.,   5.],
           [ nan,   8.,   9.],
           [  5.,   4.,   3.],
           [  2.,   1.,   1.],
           [ nan,   2.,   3.]])
    >>> # Upscale the missing_distribution to add more missing values
    >>> abs_missing_rate = 0.2
    >>> missing_distribution = {1: [0.8 * abs_missing_rate,
    ...                             0.2 * abs_missing_rate,
    ...                             0]}
    >>> vd = ValueDropper(missing_distribution=missing_distribution,
    ...                   random_state=0)
    >>> vd.transform(X, y)
    array([[  0.,   1.,   2.],
           [  3.,   4.,   5.],
           [  6.,   7.,   8.],
           [  9.,   0.,   1.],
           [  2.,   3.,   4.],
           [ nan,   9.,   8.],
           [ nan,  nan,   5.],
           [ nan,   8.,   9.],
           [ nan,   4.,   3.],
           [  2.,   1.,   1.],
           [ nan,   2.,   3.]])
    >>> # MCAR missingness
    >>> vd = ValueDropper(missing_distribution=abs_missing_rate,
    ...                   random_state=0)
    >>> vd.transform(X, y)
    array([[  0.,   1.,   2.],
           [  3.,   4.,   5.],
           [  6.,   7.,  nan],
           [  9.,   0.,   1.],
           [ nan,   3.,   4.],
           [  8.,   9.,   8.],
           [  1.,   0.,  nan],
           [  7.,   8.,   9.],
           [  5.,   4.,   3.],
           [ nan,  nan,   1.],
           [  1.,  nan,   3.]])
    >>> # Upscale the missing_distribution to add more missing values
    >>> # Explicitly set copy=False for inplace dropping of values
    >>> vd = ValueDropper(missing_distribution=2 * abs_missing_rate,
    ...                   copy=False, random_state=0)
    >>> _ = vd.transform(X, y)
    >>> X
    array([[  0.,   1.,   2.],
           [  3.,   4.,   5.],
           [ nan,   7.,  nan],
           [  9.,  nan,   1.],
           [ nan,   3.,   4.],
           [  8.,  nan,   8.],
           [  1.,   0.,  nan],
           [  7.,   8.,  nan],
           [  5.,   4.,   3.],
           [ nan,  nan,  nan],
           [ nan,  nan,   3.]])
    """

    def __init__(self, missing_values="NaN",
                 missing_distribution=None, copy=True, random_state=None):
        self.missing_values = missing_values
        self.missing_distribution = missing_distribution
        self.copy = copy
        self.random_state = random_state

    def transform(self, X, y=None):
        """Drop values from ``X`` according to the given distribution.

        Parameters
        ----------

        X : ndarray like of shape (n_features, n_samples)
            Data, in which the values must be dropped and set to
            ``missing_values``.

        y : array-like, shape = (n_samples,), optional for MCAR
            Target relative to X for classification or regression;
            When missing_distribution is not a dict (for MCAR missingness),
            ``y`` need not be passed.
        """
        # Validate missing_values and generate missing_mask
        if ((isinstance(self.missing_values, str) and
                (self.missing_values.lower() == "nan")) or
                np.isnan(self.missing_values)):
            missing_values = np.nan
        else:
            missing_values = self.missing_values

        # Don't allow pre-exising missing values in X, to simplify API
        X = check_array(X, dtype=('numeric'
                                  if isinstance(missing_values,
                                                (numbers.Integral, np.integer))
                                  else np.float),
                        copy=self.copy)

        n_samples, n_features = X.shape
        n_values = n_samples * n_features

        rng = check_random_state(self.random_state)

        # Validate y, and find type of missingness
        if isinstance(self.missing_distribution, dict):
            missing_type = 'nmar'
            # For NMAR
            # Validate and convert the missing_distribution dict into a
            # 2D probability distribution along the features and labels

            if y is None:
                raise ValueError("The missing_distribution is a dict "
                                 "but y is None. If missingness is to be "
                                 "related to the class labels, target class "
                                 "labels (y) must be passed.")

            target_type = type_of_target(y)
            if 'continuous' in target_type or 'multioutput' in target_type:
                raise ValueError("Value dropping based on the given "
                                 "distribution can be done only for single "
                                 "target which is discrete (classification "
                                 "tasks). The given target (y) is of type %s"
                                 % target_type)
            y = check_array(y, ensure_2d=False, dtype='numeric')

            le = LabelEncoder().fit(y)
            classes = le.classes_
            n_classes = classes.shape[0]

            drop_probs = np.zeros((n_classes, n_features), dtype=np.float64)

            for class_key, val in self.missing_distribution.items():
                # This will also validate incorrect values for class_key
                encoded_class_key = le.transform([class_key, ])[0]

                if isinstance(val, (np.floating, float)):
                    drop_probs[encoded_class_key, :] = (
                        val / float(n_features))
                elif isinstance(val, (np.ndarray, list, tuple)):
                    val = np.asarray(val)
                    if val.shape[0] != n_features:
                        raise ValueError("The shape of the per feature"
                                         " drop probabilities vector "
                                         "for label, %s, does not conform"
                                         " to the number of features, %d"
                                         % (class_key, n_features))
                    drop_probs[encoded_class_key, :] = val
                else:
                    raise ValueError("If missing_distribution is a dict with"
                                     " target labels as keys, the values of "
                                     "the dict should either be a single float"
                                     " or an array of shape (n_features,). "
                                     "%r was passed for class label %s"
                                     % (val, class_key))

        else:
            missing_type = 'mcar'

            # For MCAR
            # Validate and convert the missing_distribution dict into a
            # 1D probability distribution along the features

            drop_probs = np.zeros((1, n_features), dtype=np.float64)

            if isinstance(self.missing_distribution, (float, np.floating)):
                drop_probs[:] = self.missing_distribution / n_features
            elif isinstance(self.missing_distribution,
                            (list, tuple, np.ndarray)):
                missing_distribution = np.asarray(self.missing_distribution)
                if missing_distribution.shape[0] != n_features:
                    raise ValueError("The shape of the per feature "
                                     "drop probabilities vector does not "
                                     "conform to the number of features, %d"
                                     % n_features)
                drop_probs[:] = self.missing_distribution
            else:
                raise ValueError("missing_distribution must be a float or "
                                 " 1D vector (list, tuple or np.ndarray) of "
                                 "shape (n_features,) or dict of 1D vector / "
                                 "floats. %r was passed"
                                 % self.missing_distribution)

        if 1 - drop_probs.ravel().sum() <= np.finfo(float).eps:
            raise ValueError("The sum of all probabilities in the "
                             "missing_distribution should sum up to less "
                             "than 1. The sum was found to be %0.8f"
                             % drop_probs.ravel().sum())

        drop_counts = (drop_probs * n_values).astype(int)

        if missing_type == 'mcar':
            self._block_drop_missing_values(X, Ellipsis, 0, drop_counts, rng,
                                            missing_type, missing_values)
        else:
            # For classification, NMAR, consider missing distribution
            # based on class labels as subsets within data
            for i, class_i in enumerate(classes):
                self._block_drop_missing_values(X, y == class_i, i,
                                                drop_counts, rng,
                                                missing_type, missing_values)
        return X

    def _block_drop_missing_values(self, X, samples_mask, encoded_label,
                                   drop_counts, rng, missing_type,
                                   missing_values):
        """Helper to insert missing values in given block (label)"""
        n_features = X.shape[1]
        this_block_indices = np.arange(X.shape[0])[samples_mask]
        for feature in range(n_features):
            this_n_values = X[samples_mask].shape[0]
            this_required_n_missing = drop_counts[encoded_label, feature]

            if this_required_n_missing <= 0:
                continue

            if this_required_n_missing > this_n_values:
                raise ValueError("There are no more available values at "
                                 "%sfeature - %s, to drop."
                                 # For NMAR, specify the label too
                                 % ("label - %s, " % encoded_label
                                     if missing_type == 'nmar'
                                     else "", feature))

            # Shuffle and pick this_required_n_missing indices for dropping
            picked_indices = rng.permutation(
                this_block_indices)[:this_required_n_missing]

            # Drop them
            X[picked_indices, feature] = missing_values
        return X
