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

from bz2 import BZ2File
from contextlib import closing
import gzip
import io
import os.path

import numpy as np
import scipy.sparse as sp

from ._svmlight_format import _load_svmlight_file
from .. import __version__
from ..utils import atleast2d_or_csr


def load_svmlight_file(f, n_features=None, dtype=np.float64,
                       multilabel=False, zero_based="auto"):
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
    f: {str, file-like, int}
        (Path to) a file to load. If a path ends in ".gz" or ".bz2", it will
        be uncompressed on the fly. If an integer is passed, it is assumed to
        be a file descriptor. A file-like or file descriptor will not be closed
        by this function. A file-like object must be opened in binary mode.

    n_features: int or None
        The number of features to use. If None, it will be inferred. This
        argument is useful to load several files that are subsets of a
        bigger sliced dataset: each subset might not have example of
        every feature, hence the inferred shape might vary from one
        slice to another.

    multilabel: boolean, optional
        Samples may have several labels each (see
        http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multilabel.html)

    zero_based: boolean or "auto", optional
        Whether column indices in f are zero-based (True) or one-based
        (False). If set to "auto", a heuristic check is applied to determine
        this from the file contents. Both kinds of files occur "in the wild",
        but they are unfortunately not self-identifying. Using "auto" or True
        should always be safe.

    Returns
    -------
    (X, y)

    where X is a scipy.sparse matrix of shape (n_samples, n_features),
          y is a ndarray of shape (n_samples,), or, in the multilabel case,
          a list of tuples of length n_samples.

    See also
    --------
    load_svmlight_files: similar function for loading multiple files in this
    format, enforcing the same number of features/columns on all of them.
    """
    return tuple(load_svmlight_files([f], n_features, dtype, multilabel,
                                     zero_based))


def _gen_open(f):
    if isinstance(f, int):  # file descriptor
        return io.open(f, "rb", closefd=False)
    elif not isinstance(f, basestring):
        raise TypeError("expected {str, int, file-like}, got %s" % type(f))

    _, ext = os.path.splitext(f)
    if ext == ".gz":
        return gzip.open(f, "rb")
    elif ext == ".bz2":
        return BZ2File(f, "rb")
    else:
        return open(f, "rb")


def _open_and_load(f, dtype, multilabel, zero_based):
    if hasattr(f, "read"):
        return _load_svmlight_file(f, dtype, multilabel, zero_based)
    # XXX remove closing when Python 2.7+/3.1+ required
    with closing(_gen_open(f)) as f:
        return _load_svmlight_file(f, dtype, multilabel, zero_based)


def load_svmlight_files(files, n_features=None, dtype=np.float64,
                        multilabel=False, zero_based="auto"):
    """Load dataset from multiple files in SVMlight format

    This function is equivalent to mapping load_svmlight_file over a list of
    files, except that the results are concatenated into a single, flat list
    and the samples vectors are constrained to all have the same number of
    features.

    Parameters
    ----------
    files : iterable over {str, file-like, int}
        (Paths of) files to load. If a path ends in ".gz" or ".bz2", it will
        be uncompressed on the fly. If an integer is passed, it is assumed to
        be a file descriptor. File-likes and file descriptors will not be
        closed by this function. File-like objects must be opened in binary
        mode.

    n_features: int or None
        The number of features to use. If None, it will be inferred from the
        maximum column index occurring in any of the files.

    multilabel: boolean, optional
        Samples may have several labels each (see
        http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multilabel.html)

    zero_based: boolean or "auto", optional
        Whether column indices in files are zero-based (True) or one-based
        (False). If set to "auto", a heuristic check is applied to determine
        this from the files' contents. Both kinds of files occur "in the wild",
        but they are unfortunately not self-identifying. Using "auto" or True
        should always be safe.

    Returns
    -------
    [X1, y1, ..., Xn, yn]

    where each (Xi, yi) pair is the result from load_svmlight_file(files[i]).

    Rationale
    ---------
    When fitting a model to a matrix X_train and evaluating it against a
    matrix X_test, it is essential that X_train and X_test have the same
    number of features (X_train.shape[1] == X_test.shape[1]). This may not
    be the case if you load the files individually with load_svmlight_file.

    See also
    --------
    load_svmlight_file
    """
    r = [_open_and_load(f, dtype, multilabel, bool(zero_based)) for f in files]

    if zero_based is False \
     or zero_based == "auto" and all(np.min(indices) > 0
                                     for _, indices, _, _ in r):
        for _, indices, _, _ in r:
            indices -= 1

    if n_features is None:
        n_features = max(indices.max() for _, indices, _, _ in r) + 1

    result = []
    for data, indices, indptr, y in r:
        shape = (indptr.shape[0] - 1, n_features)
        result += sp.csr_matrix((data, indices, indptr), shape), y

    return result


def _dump_svmlight(X, y, f, one_based, comment):
    is_sp = int(hasattr(X, "tocsr"))
    if X.dtype == np.float64:
        value_pattern = u"%d:%0.16e"
    else:
        value_pattern = u"%d:%f"

    if y.dtype.kind == 'i':
        line_pattern = u"%d %s\n"
    else:
        line_pattern = u"%f %s\n"

    f.write("# Generated by dump_svmlight_file from scikit-learn %s\n"
            % __version__)
    f.write("# Column indices are %s-based\n" % ["zero", "one"][one_based])

    if comment:
        f.write("#\n")
        f.writelines("# %s\n" % line for line in comment)

    for i in xrange(X.shape[0]):
        s = u" ".join([value_pattern % (j + one_based, X[i, j])
                       for j in X[i].nonzero()[is_sp]])
        f.write((line_pattern % (y[i], s)).encode('ascii'))


def dump_svmlight_file(X, y, f, zero_based=True, comment=None):
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

    f : string or file-like in binary mode
        If string, specifies the path that will contain the data.
        If file-like, data will be written to f. f should be opened in binary
        mode.

    zero_based : boolean, optional
        Whether column indices should be written zero-based (True) or one-based
        (False).

    comment : string, optional
        Comment to insert at the top of the file. This should be either a
        Unicode string, which will be encoded as UTF-8, or an ASCII byte
        string.
    """
    if comment is not None:
        # Convert comment string to list of lines in UTF-8.
        # If a byte string is passed, then check whether it's ASCII;
        # if a user wants to get fancy, they'll have to decode themselves.
        # Avoid mention of str and unicode types for Python 3.x compat.
        if isinstance(comment, bytes):
            comment.decode("ascii")     # just for the exception
        else:
            comment = comment.encode("utf-8")
        if "\0" in comment:
            raise ValueError("comment string contains NUL byte")
        comment = comment.splitlines()

    y = np.asarray(y)
    if y.ndim != 1:
        raise ValueError("expected y of shape [n_samples], got %r" % y)

    X = atleast2d_or_csr(X)
    if X.shape[0] != y.shape[0]:
        raise ValueError("X.shape[0] and y.shape[0] should be the same, "
                         "got: %r and %r instead." % (X.shape[0], y.shape[0]))

    one_based = not zero_based

    if hasattr(f, "write"):
        _dump_svmlight(X, y, f, one_based, comment)
    else:
        with open(f, "wb") as f:
            _dump_svmlight(X, y, f, one_based, comment)
