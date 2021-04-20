
.. _GroupTimeSeriesSplit:

=================================================
sklearn.model_selection.GroupTimeSeriesSplit
=================================================
.. code-block:: python

   class sklearn.model_selection.GroupTimeSeriesSplit(n_splits=5, *, max_train_size=None, test_size=None, gap=0)

| *GroupTimeSeriesSplit* combines *TimeSeriesSplit* with the Group awareness of *GroupKFold*.
|
| Like *TimeSeriesSplit* this  also returns first *k* folds as train set and the *(k+1)* th fold as test set.
|
| Since the Group applies on this class, the same group will not appear in two different
 folds(the number of distinct groups has to be at least equal to the number of folds) which make sure the i.i.d. assumption will not be broken.

| All operations of this CV strategy are done at the group level.
| So all our parameters, not limited to splits, including test_size, gap, and max_train_size, all represent the constraints on the number of groups.


Parameters: 
-----------
| **n_splits;int,default=5**
|
|   Number of splits. Must be at least 2.
|
| **max_train_size:int, default=None**
|
|   Maximum number of groups for a single training set.
|
| **test_size:int, default=None**
|
|   Used to limit the number of groups in the test set. Defaults to ``n_samples // (n_splits + 1)``, which is the maximum allowed value with ``gap=0``.
|
| **gap:int, default=0**
|
|  Number of groups in samples to exclude from the end of each train set before the test set.

Example 1:
---------
.. code-block:: python

>>> import numpy as np
>>> from sklearn.model_selection import GroupTimeSeriesSplit
>>> groups = np.array(['a', 'a', 'a', 'a', 'a', 'a',
...                    'b', 'b', 'b', 'b', 'b',
...                    'c', 'c', 'c', 'c',
...                    'd', 'd', 'd'])
>>> gtss = GroupTimeSeriesSplit(n_splits=3)
>>> for train_idx, test_idx in gtss.split(groups, groups=groups):
...     print("TRAIN:", train_idx, "TEST:", test_idx)
...     print("TRAIN GROUP:", groups[train_idx],
...           "TEST GROUP:", groups[test_idx])
TRAIN: [0, 1, 2, 3, 4, 5] TEST: [6, 7, 8, 9, 10]
TRAIN GROUP: ['a' 'a' 'a' 'a' 'a' 'a']
TEST GROUP: ['b' 'b' 'b' 'b' 'b']
TRAIN: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] TEST: [11, 12, 13, 14]
TRAIN GROUP: ['a' 'a' 'a' 'a' 'a' 'a' 'b' 'b' 'b' 'b' 'b']
TEST GROUP: ['c' 'c' 'c' 'c']
TRAIN: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
TEST: [15, 16, 17]
TRAIN GROUP: ['a' 'a' 'a' 'a' 'a' 'a' 'b' 'b' 'b' 'b' 'b' 'c' 'c' 'c' 'c']
TEST GROUP: ['d' 'd' 'd']

Example 2:
---------
.. code-block:: python

>>> import numpy as np
>>> from sklearn.model_selection import GroupTimeSeriesSplit
>>> groups = np.array(['a', 'a', 'a', 'a', 'a', 'a',\
                       'b', 'b', 'b', 'b', 'b',\
                       'c', 'c', 'c', 'c',\
                       'd', 'd', 'd'])
>>> gtss = GroupTimeSeriesSplit(n_splits=2, test_size=1, gap=1,\
                                max_train_size=3)
>>> for train_idx, test_idx in gtss.split(groups, groups=groups):
...     print("TRAIN:", train_idx, "TEST:", test_idx)
...     print("TRAIN GROUP:", groups[train_idx],\
              "TEST GROUP:", groups[test_idx])
TRAIN: [0, 1, 2, 3, 4, 5] TEST: [11, 12, 13, 14]
TRAIN GROUP: ['a' 'a' 'a' 'a' 'a' 'a'] TEST GROUP: ['c' 'c' 'c' 'c']
TRAIN: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10] TEST: [15, 16, 17]
TRAIN GROUP: ['a' 'a' 'a' 'a' 'a' 'a' 'b' 'b' 'b' 'b' 'b']
TEST GROUP: ['d' 'd' 'd']

Methods: 
--------
| **get_n_splits([X, y, groups])**
|
|   Returns the number of splitting iterations in the cross-validator
|   *Parameters:*
|       *X: object*
|           Always ignored, exists for compatibility.
|       *y: object*
|           Always ignored, exists for compatibility.
|       *groups: object*
|           Always ignored, exists for compatibility.
|   *Returns:*
|       *n_splits: int*
|           Returns the number of splitting iterations in the cross-validator.
|
| **split(X[groups, y])**
|
|   Generate indices to split data into training and test set by group.
|   *Parameters:*
|       *X : array-like of shape (n_samples, n_features)*
|            Training data, where n_samples is the number of samples
|            and n_features is the number of features.
|       *y : array-like of shape (n_samples,)*
|            Always ignored, exists for compatibility.
|       *groups : array-like of shape (n_samples,)*
|            Group labels for the samples used while splitting the dataset into
|            train/test set.
|   *Yields:*
|       *train : ndarray*
|            The training set indices for that split.
|       *test : ndarray*
|            The testing set indices for that split.

