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

    If the missing type is NMAR, a ``missing_proba`` parameter can be passed
    to drop values conforming to the given drop-probabilities.


    Parameters
    ----------

    missing_values : {"NaN" (or np.nan) | int}, default "NaN"
        The value to insert to indicate missingness.

    missing_proba : dict of floats or dict of vector of floats
        If ``missing_proba`` is a float within range [0, 1), it represents the
        probability with which the values will be dropped.

        The values are dropped (approximately) uniformly across all labels and
        features. This type of missingness is referred to as MCAR.

        To vary the proportion of values dropped across each feature,
        individual drop-probabilities for each feature can be specified as a 1D
        vector of shape ``(n_features,)``.

        If missingness is not MCAR, ``missing_proba`` can be used to specify
        the drop-probabilities on a per-label (and if needed further on
        per-feature basis.).

        If ``missing_proba`` is a dict of floats::
            {1: 0.2, 2: 0.3, 3: 0.5}

        This represents, the drop-probabilities for samples of each
        class-label. The missing values are evenly spread across all the
        features.

        If ``missing_proba`` is a dict of vectors (and scalars)::

            {0: 0.1,
             3: [0.1, 0.15]}

        Note that the shape of the vector must be ``(n_features,)``

        Samples from class 0 are dropped with probability of 0.1 for each
        feature and those from class 3 are dropped with a probability of 0.1
        in feature 0, 0.15 in feature 1 while there are no values dropped from
        samples of class 1 and 2.

        Note that the samples are randomly chosen "*for each feature*".

    copy : bool, default False
        Whether to copy the data or work inplace.

    random_state : int, optional
        The seed for the numpy's random number generator.

        If ``random_state`` is set to an integer, the ``missing_proba``
        can be scaled uniformly with the assumption that all the values
        dropped with a smaller scale will exist in the larger scaled version::
            missing_proba_25pc = {0: 0.05, 3: [0.05, 0.075, 0.075]}
            missing_proba_50pc = {0: 0.1, 3: [0.1, 0.15, 0.15]}

        The missing values dropped with ``missing_proba_25pc`` will also
        exist in ``missing_proba_50pc``.

        This guarantee does not apply when relative probabilities
        within the ``missing_proba`` change between two settings::
            missing_proba_25pc = {0: 0.05, 3: [0.05, 0.075, 0.075]}
            # The below dist. is not a scaled version of above
            missing_proba_50pc = {0: 0.09, 3: [0.11, 0.15, 0.15]}


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
    >>> missing_proba = {1: [0.8 * abs_missing_rate,
    ...                      0.2 * abs_missing_rate,
    ...                      0]}
    >>> vd = ValueDropper(missing_proba=missing_proba, random_state=0)
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
    >>> # Upscale the missing_proba to add more missing values
    >>> abs_missing_rate = 0.2
    >>> missing_proba = {1: [0.8 * abs_missing_rate,
    ...                             0.2 * abs_missing_rate,
    ...                             0]}
    >>> vd = ValueDropper(missing_proba=missing_proba, random_state=0)
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
    >>> vd = ValueDropper(missing_proba=abs_missing_rate, random_state=0)
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
    >>> # Upscale the missing_proba to add more missing values
    >>> # Explicitly set copy=False for inplace dropping of values
    >>> vd = ValueDropper(missing_proba=2 * abs_missing_rate,
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
                 missing_proba=None, copy=True, random_state=None):
        self.missing_values = missing_values
        self.missing_proba = missing_proba
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
            When missing_proba is not a dict (for MCAR missingness),
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
        rng = check_random_state(self.random_state)

        # Validate y, and find type of missingness
        if isinstance(self.missing_proba, dict):
            missing_type = 'nmar'
            # For NMAR
            # Validate and convert the missing_proba dict into a
            # 2D probability distribution along the features and labels

            if y is None:
                raise ValueError("The missing_proba is a dict "
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

            for class_key, val in self.missing_proba.items():
                # This will also validate incorrect values for class_key
                encoded_class_key = le.transform([class_key, ])[0]

                if isinstance(val, (np.ndarray, list, tuple)):
                    val = np.asarray(val)
                    if val.shape[0] != n_features:
                        raise ValueError("The shape of the per feature "
                                         "drop-probabilities vector "
                                         "for label, %s, does not conform "
                                         "to the number of features, %d"
                                         % (class_key, n_features))
                elif not isinstance(val, (np.floating, float,
                                          numbers.Integral, np.integer)):
                    raise ValueError("If missing_proba is a dict with "
                                     "target labels as keys, the values of "
                                     "the dict should either be a single "
                                     "float or an array of shape "
                                     "(n_features,). %r was passed for class "
                                     "label %s" % (val, class_key))

                drop_probs[encoded_class_key, :] = val

        else:
            missing_type = 'mcar'
            # For MCAR
            # Validate and convert the missing_proba dict into a
            # 1D probability distribution along the features

            drop_probs = np.zeros((1, n_features), dtype=np.float64)

            if isinstance(self.missing_proba, (list, tuple, np.ndarray)):
                # Convert to ndarray and check shape
                missing_proba = np.asarray(self.missing_proba)
                if missing_proba.shape[0] != n_features:
                    raise ValueError("The shape of the per feature "
                                     "drop-probabilities vector does not "
                                     "conform to the number of features, %d"
                                     % n_features)
            elif not isinstance(self.missing_proba,
                                (np.floating, float, numbers.Integral,
                                 np.integer)):
                raise ValueError("missing_proba must be a float or "
                                 "1D vector (list, tuple or np.ndarray) of "
                                 "shape (n_features,) or dict of 1D vector / "
                                 "floats. %r was passed"
                                 % self.missing_proba)

            drop_probs[:] = self.missing_proba
            # Hack to simplify code
            classes = [0, ]
            y = np.zeros(n_samples)

        if np.any(drop_probs < 0) or np.any(drop_probs > 1):
            raise ValueError("All the individual drop-probabilities should be "
                             "within the range of [0, 1]. The given "
                             "missing_proba does not conform to that. %r"
                             % self.missing_proba)

        for i, class_i in enumerate(classes):
            samples_mask = (y == class_i)
            this_n_samples = samples_mask.sum()
            this_block_indices = np.arange(n_samples)[samples_mask]

            for feature in range(n_features):
                # Shuffle even if this_required_n_missing is 0, to maintain
                # consistency in generated missing values for successively
                # increasing % of missing values.%
                shuffled_indices = rng.permutation(this_block_indices)
                this_required_n_missing = int(round(drop_probs[i, feature] *
                                                    this_n_samples))
                if this_required_n_missing == 0:
                    continue

                if this_required_n_missing > this_n_samples:
                    raise ValueError("There are no more available values at "
                                     "%sfeature - %s, to drop."
                                     # For NMAR, specify the label too
                                     % ("class label - %s, " % class_i
                                        if missing_type == 'nmar' else "",
                                        feature))

                # Drop them
                X[shuffled_indices[:this_required_n_missing],
                  feature] = missing_values

        return X
