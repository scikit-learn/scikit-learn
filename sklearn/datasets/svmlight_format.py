""" This module implements a loader for the svmlight / libsvm
sparse dataset format.  """

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
# License: Simple BSD.

import numpy as np
import scipy.sparse as sp
import re

def _load_svmlight_file(file_path, buffer_mb):
    data = []
    indptr = []
    indices = []
    labels = []

    pattern = re.compile(r'\s+')

    for line in open(file_path):
        line = line.strip()

        if line.startswith("#"):
            continue

        # Remove inline comments.
        line = line.split("#")
        line = line[0].strip()

        line = re.sub(pattern, ' ', line)
        y, features = line.split(" ", 1)

        labels.append(float(y))
        indptr.append(len(data))

        for feat in features.split(" "):
            idx, value = feat.split(":")
            indices.append(int(idx))
            data.append(float(value))

    indptr.append(len(data))

    return np.array(data, dtype=np.double), \
           np.array(indices, dtype=np.int), \
           np.array(indptr, dtype=np.int), \
           np.array(labels, dtype=np.double)


def load_svmlight_file(file_path, other_file_path=None,
                         n_features=None, buffer_mb=40):
    """Load datasets in the svmlight / libsvm format directly into
    scipy sparse CSR matrices.

    Parameters
    ----------
    file_path: str
        Path to a file to load.

    other_file_path: str or None
        Path to another file to load. scikit-learn will make sure that the
        number of features in the returned matrix is the same as for
        file_path.

    n_features: int or None
        The number of features to use. If None, it will be inferred.

    buffer_mb: int (default: 40)
        The size of the buffer used while loading the dataset in mega-bytes.

    Returns
    -------
        (X, y)

        where X is a scipy.sparse matrix of shape (n_samples, n_features),
              y is a ndarray of shape (n_samples,),

        or, if other_file_path is not None,

        (X1, y1, X2, y2)

        where X1 and X2 are scipy.sparse matrices of shape
                            (n_samples1, n_features) and
                            (n_samples2, n_features),
              y1 and y2 are ndarrays of shape (n_samples1,) and (n_samples2,).

    Note
    ----
    When fitting a model to a matrix X_train and evaluating it against a matrix
    X_test, it is essential that X_train and X_test have the same number of
    features (X_train.shape[1] == X_test.shape[1]). This may not be the case if
    you load them with load_svmlight_format separately. To address this
    problem, we recommend to use load_svmlight_format(train_file, test_file)
    or load_svmlight_format(test_file, n_features=X_train.shape[1]).
    """
    data, indices, indptr, labels = _load_svmlight_file(file_path, buffer_mb)

    if n_features is not None:
        shape = (indptr.shape[0] - 1, n_features)
    else:
        shape = None    # inferred

    X_train = sp.csr_matrix((data, indices, indptr), shape)

    ret = [X_train, labels]

    if other_file_path is not None:
        tup = _load_svmlight_file(other_file_path, buffer_mb)
        data, indices, indptr, labels = tup

        if n_features is None:
            n_features = X_train.shape[1]

        shape = (indptr.shape[0] - 1, n_features)

        X_test = sp.csr_matrix((data, indices, indptr), shape)

        ret.append(X_test)
        ret.append(labels)

    return tuple(ret)
