"""This module implements a loader and dumper for the svmlight format

This format is a text-based format, with one sample per line. It does
not store zero valued features hence is suitable for sparse dataset.

The first element of each line can be used to store a target variable to
predict.

This format is used as the default format for both svmlight and the
libsvm command line programs.
"""

# Authors: Mathieu Blondel <mathieu@mblondel.org>
#          Lars Buitinck <L.J.Buitinck@uva.nl>
# License: Simple BSD.

import numpy as np
import scipy.sparse as sp
import re

def _load_svmlight_file(f, buffer_mb, n_features):
    data = []
    indptr = []
    indices = []
    labels = []

    pattern = re.compile(r'\s+')

    for line in f:
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
    indptr = np.array(indptr, dtype=np.int)

    if n_features is not None:
        shape = (indptr.shape[0] - 1, n_features)
    else:
        shape = None    # inferred

    X = sp.csr_matrix((np.array(data, dtype=np.double),
                       np.array(indices, dtype=np.int),
                       indptr),
                      shape)

    return X, np.array(labels, dtype=np.double)


def _load_svmlight(f, *args):
    if hasattr(f, "read"):
        return _load_svmlight_file(f, *args)
    with open(f) as f:
        return _load_svmlight_file(f, *args)


def load_svmlight_file(f, other_file=None, n_features=None, buffer_mb=40):
    """Load datasets in the svmlight / libsvm format into sparse CSR matrix

    This format is a text-based format, with one sample per line. It does
    not store zero valued features hence is suitable for sparse dataset.

    The first element of each line can be used to store a target variable
    to predict.

    This format is used as the default format for both svmlight and the
    libsvm command line programs.

    Parsing a text based source can be expensive. When working on
    repeatedly on the same dataset, it is recommended to wrap this
    loader with joblib.Memory.cache to store a memmapped backup of the
    CSR results of the first call and benefit from the near instantaneous
    loading of memmapped structures for the subsequent calls.

    This implementation is naive: it does allocate too much memory and
    is slow since written in python. On large datasets it is recommended
    to use an optimized loader such as:

      https://github.com/mblondel/svmlight-loader

    Parameters
    ----------
    f: str or file-like
        (Path to) a file to load.

    other_file: str or file-like, optional
        (Path to) another file to load. The benefit over calling this function
        twice for the files is that n_features is enforced on both datasets.

    n_features: int or None
        The number of features to use. If None, it will be inferred. This
        argument is useful to load several files that are subsets of a
        bigger sliced dataset: each subset might not have example of
        every feature, hence the inferred shape might vary from one
        slice to another.

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
    ret = _load_svmlight(f, buffer_mb, n_features)

    if other_file is not None:
        ret += _load_svmlight(other_file, buffer_mb, n_features)

    return ret


def _dump_svmlight(X, y, f):
    if X.shape[0] != y.shape[0]:
        raise ValueError("X.shape[0] and y.shape[0] should be the same.")

    is_sp = int(hasattr(X, "tocsr"))

    for i in xrange(X.shape[0]):
        s = " ".join(["%d:%f" % (j, X[i, j]) for j in X[i].nonzero()[is_sp]])
        f.write("%f %s\n" % (y[i], s))


def dump_svmlight_file(X, y, f):
    """Dump the dataset in svmlight / libsvm file format.

    This format is a text-based format, with one sample per line. It does
    not store zero valued features hence is suitable for sparse dataset.

    The first element of each line can be used to store a target variable
    to predict.

    Parameters
    ----------
    X : {array-like, sparse matrix}, shape = [n_samples, n_features]
        Training vectors, where n_samples is the number of samples and
        n_features is the number of features.

    y : array-like, shape = [n_samples]
        Target values.

    f : str or file-like
    """
    if hasattr(f, "write"):
        _dump_svmlight(X, y, f)
    else:
        with open(f, "w") as f:
            _dump_svmlight(X, y, f)
