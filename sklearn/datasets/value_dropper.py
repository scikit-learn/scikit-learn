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

    missing_values : {"NaN" (or np.nan) | int | float}, default "NaN"
        The value to insert to indicate missingness.

    missing_proba : dict of floats or dict of vector of floats
        To vary the proportion of values dropped across each feature,
        individual drop-probabilities for each feature can be specified as a 1D
        array-like of shape (n_features, ) (e.g. [0.1, 0.15, 0.1]).

        If missingness is not MCAR, a dict of floats can be used to specify
        the drop-probabilities on a per-label basis
        (e.g. {1: 0.2, 2: 0.3, 3: 0.5}).

        This dict can also contains some 1D array-like of shape (n_features, )
        to vary drop-probabilities across features
        (e.g. {1: 0.1, 3: [0.1, 0.15, 0.1]}).

    copy : bool, default False
        Whether to copy the data or work inplace.

    random_state : int, optional
        The seed for the numpy's random number generator.

        If ``random_state`` is set to an integer, the ``missing_proba``
        can be upscaled safely with the assumption that all the values
        dropped with a smaller scale will exist in the larger scaled version::
            missing_proba_1 = {0: 0.1, 3: [0.3, 0.1, 0.1]}
            missing_proba_2 = {0: 0.1, 1:0.2, 3: [0.6, 0.1, 0.8]}

        The missing values dropped with ``missing_proba_1`` will also
        be dropped with ``missing_proba_2``.


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
    >>> # NMAR missingness -
    >>> # Drop values from samples of class 1 alone based on the below
    >>> # missing_proba hence making it Not Missing At Random missingness.
    >>> missing_proba = {1: [0.2,  # Drop 20% values from feature 0 for class 0
    ...                      0.2,  # and class 1
    ...                      0]}   # Do not drop any values from feature 2
    >>> vd = ValueDropper(missing_proba=missing_proba, random_state=0)
    >>> vd.transform(X, y)
    array([[  0.,   1.,   2.],
           [  3.,   4.,   5.],
           [  6.,   7.,   8.],
           [  9.,   0.,   1.],
           [  2.,   3.,   4.],
           [  8.,   9.,   8.],
           [  1.,   0.,   5.],
           [  7.,   8.,   9.],
           [ nan,   4.,   3.],
           [  2.,   1.,   1.],
           [  1.,  nan,   3.]])
    >>> # Upscale the missing_proba to add more missing values in feature 0
    >>> # Also add a few missing values in all features for class 0 samples.
    >>> missing_proba = {1: [0.4, 0.2, 0], 0: 0.6}
    >>> vd = ValueDropper(missing_proba=missing_proba, random_state=0)
    >>> vd.transform(X, y)
    array([[ nan,  nan,   2.],
           [ nan,  nan,  nan],
           [ nan,  nan,   8.],
           [  9.,   0.,  nan],
           [  2.,   3.,  nan],
           [  8.,   9.,   8.],
           [  1.,   0.,   5.],
           [  7.,   8.,   9.],
           [ nan,   4.,   3.],
           [  2.,   1.,   1.],
           [ nan,  nan,   3.]])
    >>> # MCAR missingness -
    >>> # 30% of values in each feature Missing Completely At Random
    >>> vd = ValueDropper(missing_proba=0.3, random_state=0)
    >>> vd.transform(X, y)
    array([[  0.,   1.,   2.],
           [  3.,   4.,   5.],
           [ nan,   7.,  nan],
           [  9.,  nan,   1.],
           [ nan,   3.,   4.],
           [  8.,   9.,   8.],
           [  1.,   0.,  nan],
           [  7.,   8.,  nan],
           [  5.,   4.,   3.],
           [ nan,  nan,   1.],
           [  1.,  nan,   3.]])
    >>> # Upscale the missing_proba to add more missing values in feature 0 and
    >>> # 1 alone. Retain the same drop-probability for feature 2
    >>> # Explicitly set copy=False for inplace dropping of values
    >>> vd = ValueDropper(missing_proba=[0.6, 0.8, 0.3],
    ...                   copy=False, random_state=0)
    >>> _ = vd.transform(X, y)
    >>> X
    array([[  0.,  nan,   2.],
           [ nan,  nan,   5.],
           [ nan,  nan,  nan],
           [  9.,  nan,   1.],
           [ nan,  nan,   4.],
           [  8.,  nan,   8.],
           [ nan,   0.,  nan],
           [ nan,   8.,  nan],
           [  5.,  nan,   3.],
           [ nan,  nan,   1.],
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

        X : array-like of shape (n_features, n_samples)
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
            # For NMAR
            # Validate and convert the missing_proba dict into a
            # 2D probability distribution along the features and labels
            missing_type = 'nmar'

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

            class_keys, probas = zip(*self.missing_proba.items())
            encoded_class_keys = le.transform(class_keys)
        else:
            # For MCAR
            # Validate and convert the missing_proba dict into a
            # 1D probability distribution along the features
            missing_type = 'mcar'

            drop_probs = np.zeros((1, n_features), dtype=np.float64)

            # Hack to simplify and unify missing generation code for nmar/mcar
            classes = class_keys = encoded_class_keys = (0, )
            probas = (self.missing_proba, )
            y = np.zeros(n_samples)

        # For both nmar/mcar
        for encoded_class_key, class_key, proba in zip(encoded_class_keys,
                                                       class_keys, probas):
            if isinstance(proba, (np.ndarray, list, tuple)):
                proba = np.asarray(proba)
                if proba.shape[0] != n_features:
                    raise ValueError("%s shape of the per feature "
                                     "drop-probabilities vector "
                                     "does not conform to the number of "
                                     "features, %d"
                                     % ("For label, %s, the" % class_key
                                        if missing_type == 'nmar'
                                        else "The", n_features))
            elif not isinstance(proba, (np.floating, float,
                                      numbers.Integral, np.integer)):
                raise ValueError("%s value must be a float or "
                                 "1D vector (list, tuple or np.ndarray) of "
                                 "shape (n_features,)%s %r was passed."
                                 % ("For label, %s, probability" % class_key
                                    if missing_type == 'nmar'
                                    else 'Probability',
                                    " or dict of floats/1D vectors."
                                    if missing_type == 'mcar' else "", proba))

            drop_probs[encoded_class_key, :] = proba

        if np.any(drop_probs < 0) or np.any(drop_probs > 1):
            raise ValueError("All the individual drop-probabilities should be "
                             "within the range of [0, 1]. The given "
                             "missing_proba does not conform to that. %r"
                             % self.missing_proba)

        # Generate random_states for each feature / label in advance
        # This is important to maintain consistency in generated missing values
        # for successively increasing missing percent.
        random_states = rng.randint(0, np.iinfo(np.int32).max,
                                    drop_probs.shape)

        for i, class_i in enumerate(classes):
            samples_mask = (y == class_i)
            this_n_samples = samples_mask.sum()
            this_block_indices = np.arange(n_samples)[samples_mask]

            for feature in range(n_features):
                this_required_n_missing = int(round(drop_probs[i, feature] *
                                                    this_n_samples))
                if this_required_n_missing == 0:
                    continue

                this_rng = check_random_state(random_states[i, feature])
                shuffled_indices = this_rng.permutation(this_block_indices)

                # Drop them
                X[shuffled_indices[:this_required_n_missing],
                  feature] = missing_values

        return X
