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
#          Olivier Grisel <olivier.grisel@ensta.org>
# License: Simple BSD.

import numpy as np
from ._svmlight_format import _load_svmlight_file


def load_svmlight_file(f, n_features=None, dtype=np.float64,
                       multilabel=False):
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
    f: str or file-like open in binary mode.
        (Path to) a file to load.

    n_features: int or None
        The number of features to use. If None, it will be inferred. This
        argument is useful to load several files that are subsets of a
        bigger sliced dataset: each subset might not have example of
        every feature, hence the inferred shape might vary from one
        slice to another.

    multilabel: boolean, optional
        Samples may have several labels each (see
        http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multilabel.html)

    Returns
    -------
    (X, y)

    where X is a scipy.sparse matrix of shape (n_samples, n_features),
          y is a ndarray of shape (n_samples,), or, in the multilabel case,
          a list of tuples of length n_samples.
    """
    if hasattr(f, "read"):
        return _load_svmlight_file(f, n_features, dtype, multilabel)
    with open(f, 'rb') as f:
        return _load_svmlight_file(f, n_features, dtype, multilabel)


def load_svmlight_files(files, n_features=None, dtype=np.float64,
                        multilabel=False):
    """Load dataset from multiple files in SVMlight format

    This function is equivalent to mapping load_svmlight_file over a list of
    files, except that the results are concatenated into a single, flat list
    and the samples vectors are constrained to all have the same number of
    features.

    Parameters
    ----------
    files : iterable over {str, file-like}
        (Paths to) files to load.

    n_features: int or None
        The number of features to use. If None, it will be inferred from the
        first file. This argument is useful to load several files that are
        subsets of a bigger sliced dataset: each subset might not have
        examples of every feature, hence the inferred shape might vary from
        one slice to another.

    multilabel: boolean, optional
        Samples may have several labels each (see
        http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multilabel.html)

    Returns
    -------
    [X1, y1, ..., Xn, yn]

    where each (Xi, yi) pair is the result from load_svmlight_file(files[i]).

    Rationale
    ---------
    When fitting a model to a matrix X_train and evaluating it against a
    matrix X_test, it is essential that X_train and X_test have the same
    number of features (X_train.shape[1] == X_test.shape[1]). This may not
    be the case if you load them with load_svmlight_file separately.

    See also
    --------
    load_svmlight_file
    """
    files = iter(files)
    result = list(load_svmlight_file(files.next(), n_features, dtype))
    n_features = result[0].shape[1]

    for f in files:
        result += load_svmlight_file(f, n_features, dtype, multilabel)

    return result


def _dump_svmlight(X, y, f):
    if X.shape[0] != y.shape[0]:
        raise ValueError("X.shape[0] and y.shape[0] should be the same, "
                         "got: %r and %r instead." % (X.shape[0], y.shape[0]))

    is_sp = int(hasattr(X, "tocsr"))

    for i in xrange(X.shape[0]):
        s = u" ".join([u"%d:%f" % (j + 1, X[i, j])
                       for j in X[i].nonzero()[is_sp]])
        f.write((u"%f %s\n" % (y[i], s)).encode('ascii'))


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

    f : str or file-like in binary mode
    """
    if hasattr(f, "write"):
        _dump_svmlight(X, y, f)
    else:
        with open(f, "wb") as f:
            _dump_svmlight(X, y, f)
