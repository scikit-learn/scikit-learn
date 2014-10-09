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
# License: BSD 3 clause

from contextlib import closing
import io
import os.path

import numpy as np
import scipy.sparse as sp

from ._svmlight_format import _load_svmlight_file
from .. import __version__
from ..externals import six
from ..externals.six import u, b
from ..externals.six.moves import range, zip
from ..utils import check_array
from ..utils.fixes import frombuffer_empty


def load_svmlight_file(f, n_features=None, dtype=np.float64,
                       multilabel=False, zero_based="auto", query_id=False):
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

    In case the file contains a pairwise preference constraint (known
    as "qid" in the svmlight format) these are ignored unless the
    query_id parameter is set to True. These pairwise preference
    constraints can be used to constraint the combination of samples
    when using pairwise loss functions (as is the case in some
    learning to rank problems) so that only pairs with the same
    query_id value are considered.

    This implementation is written in Cython and is reasonably fast.
    However, a faster API-compatible loader is also available at:

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
        bigger sliced dataset: each subset might not have examples of
        every feature, hence the inferred shape might vary from one
        slice to another.

    multilabel: boolean, optional
        Samples may have several labels each (see
        http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multilabel.html)

    zero_based: boolean or "auto", optional
        Whether column indices in f are zero-based (True) or one-based
        (False). If column indices are one-based, they are transformed to
        zero-based to match Python/NumPy conventions.
        If set to "auto", a heuristic check is applied to determine this from
        the file contents. Both kinds of files occur "in the wild", but they
        are unfortunately not self-identifying. Using "auto" or True should
        always be safe.

    query_id: boolean, defaults to False
        If True, will return the query_id array for each file.

    Returns
    -------
    X: scipy.sparse matrix of shape (n_samples, n_features)

    y: ndarray of shape (n_samples,), or, in the multilabel a list of
        tuples of length n_samples.

    query_id: array of shape (n_samples,)
       query_id for each sample. Only returned when query_id is set to
       True.

    See also
    --------
    load_svmlight_files: similar function for loading multiple files in this
    format, enforcing the same number of features/columns on all of them.
    """
    return tuple(load_svmlight_files([f], n_features, dtype, multilabel,
                                     zero_based, query_id))


def _gen_open(f):
    if isinstance(f, int):  # file descriptor
        return io.open(f, "rb", closefd=False)
    elif not isinstance(f, six.string_types):
        raise TypeError("expected {str, int, file-like}, got %s" % type(f))

    _, ext = os.path.splitext(f)
    if ext == ".gz":
        import gzip
        return gzip.open(f, "rb")
    elif ext == ".bz2":
        from bz2 import BZ2File
        return BZ2File(f, "rb")
    else:
        return open(f, "rb")


def _open_and_load(f, dtype, multilabel, zero_based, query_id):
    if hasattr(f, "read"):
        actual_dtype, data, ind, indptr, labels, query = \
            _load_svmlight_file(f, dtype, multilabel, zero_based, query_id)
    # XXX remove closing when Python 2.7+/3.1+ required
    else:
        with closing(_gen_open(f)) as f:
            actual_dtype, data, ind, indptr, labels, query = \
                _load_svmlight_file(f, dtype, multilabel, zero_based, query_id)

    # convert from array.array, give data the right dtype
    if not multilabel:
        labels = frombuffer_empty(labels, np.float64)
    data = frombuffer_empty(data, actual_dtype)
    indices = frombuffer_empty(ind, np.intc)
    indptr = np.frombuffer(indptr, dtype=np.intc)   # never empty
    query = frombuffer_empty(query, np.intc)

    data = np.asarray(data, dtype=dtype)    # no-op for float{32,64}
    return data, indices, indptr, labels, query


def load_svmlight_files(files, n_features=None, dtype=np.float64,
                        multilabel=False, zero_based="auto", query_id=False):
    """Load dataset from multiple files in SVMlight format

    This function is equivalent to mapping load_svmlight_file over a list of
    files, except that the results are concatenated into a single, flat list
    and the samples vectors are constrained to all have the same number of
    features.

    In case the file contains a pairwise preference constraint (known
    as "qid" in the svmlight format) these are ignored unless the
    query_id parameter is set to True. These pairwise preference
    constraints can be used to constraint the combination of samples
    when using pairwise loss functions (as is the case in some
    learning to rank problems) so that only pairs with the same
    query_id value are considered.

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

        This can be set to a higher value than the actual number of features
        in any of the input files, but setting it to a lower value will cause
        an exception to be raised.

    multilabel: boolean, optional
        Samples may have several labels each (see
        http://www.csie.ntu.edu.tw/~cjlin/libsvmtools/datasets/multilabel.html)

    zero_based: boolean or "auto", optional
        Whether column indices in f are zero-based (True) or one-based
        (False). If column indices are one-based, they are transformed to
        zero-based to match Python/NumPy conventions.
        If set to "auto", a heuristic check is applied to determine this from
        the file contents. Both kinds of files occur "in the wild", but they
        are unfortunately not self-identifying. Using "auto" or True should
        always be safe.

    query_id: boolean, defaults to False
        If True, will return the query_id array for each file.

    Returns
    -------
    [X1, y1, ..., Xn, yn]
    where each (Xi, yi) pair is the result from load_svmlight_file(files[i]).

    If query_id is set to True, this will return instead [X1, y1, q1,
    ..., Xn, yn, qn] where (Xi, yi, qi) is the result from
    load_svmlight_file(files[i])

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
    r = [_open_and_load(f, dtype, multilabel, bool(zero_based), bool(query_id))
         for f in files]

    if (zero_based is False
            or zero_based == "auto" and all(np.min(tmp[1]) > 0 for tmp in r)):
        for ind in r:
            indices = ind[1]
            indices -= 1

    n_f = max(ind[1].max() for ind in r) + 1
    if n_features is None:
        n_features = n_f
    elif n_features < n_f:
        raise ValueError("n_features was set to {},"
                         " but input file contains {} features"
                         .format(n_features, n_f))

    result = []
    for data, indices, indptr, y, query_values in r:
        shape = (indptr.shape[0] - 1, n_features)
        X = sp.csr_matrix((data, indices, indptr), shape)
        X.sort_indices()
        result += X, y
        if query_id:
            result.append(query_values)

    return result


def _dump_svmlight(X, y, f, one_based, comment, query_id):
    is_sp = int(hasattr(X, "tocsr"))
    if X.dtype.kind == 'i':
        value_pattern = u("%d:%d")
    else:
        value_pattern = u("%d:%.16g")

    if y.dtype.kind == 'i':
        line_pattern = u("%d")
    else:
        line_pattern = u("%.16g")

    if query_id is not None:
        line_pattern += u(" qid:%d")
    line_pattern += u(" %s\n")

    if comment:
        f.write(b("# Generated by dump_svmlight_file from scikit-learn %s\n"
                % __version__))
        f.write(b("# Column indices are %s-based\n"
                  % ["zero", "one"][one_based]))

        f.write(b("#\n"))
        f.writelines(b("# %s\n" % line) for line in comment.splitlines())

    for i in range(X.shape[0]):
        if is_sp:
            span = slice(X.indptr[i], X.indptr[i + 1])
            row = zip(X.indices[span], X.data[span])
        else:
            nz = X[i] != 0
            row = zip(np.where(nz)[0], X[i, nz])

        s = " ".join(value_pattern % (j + one_based, x) for j, x in row)
        if query_id is not None:
            feat = (y[i], query_id[i], s)
        else:
            feat = (y[i], s)
        f.write((line_pattern % feat).encode('ascii'))


def dump_svmlight_file(X, y, f, zero_based=True, comment=None, query_id=None):
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
        If a comment is given, then it will be preceded by one that identifies
        the file as having been dumped by scikit-learn. Note that not all
        tools grok comments in SVMlight files.

    query_id : array-like, shape = [n_samples]
        Array containing pairwise preference constraints (qid in svmlight
        format).
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
        if six.b("\0") in comment:
            raise ValueError("comment string contains NUL byte")

    y = np.asarray(y)
    if y.ndim != 1:
        raise ValueError("expected y of shape (n_samples,), got %r"
                         % (y.shape,))

    Xval = check_array(X, accept_sparse='csr')
    if Xval.shape[0] != y.shape[0]:
        raise ValueError("X.shape[0] and y.shape[0] should be the same, got"
                         " %r and %r instead." % (Xval.shape[0], y.shape[0]))

    # We had some issues with CSR matrices with unsorted indices (e.g. #1501),
    # so sort them here, but first make sure we don't modify the user's X.
    # TODO We can do this cheaper; sorted_indices copies the whole matrix.
    if Xval is X and hasattr(Xval, "sorted_indices"):
        X = Xval.sorted_indices()
    else:
        X = Xval
        if hasattr(X, "sort_indices"):
            X.sort_indices()

    if query_id is not None:
        query_id = np.asarray(query_id)
        if query_id.shape[0] != y.shape[0]:
            raise ValueError("expected query_id of shape (n_samples,), got %r"
                             % (query_id.shape,))

    one_based = not zero_based

    if hasattr(f, "write"):
        _dump_svmlight(X, y, f, one_based, comment, query_id)
    else:
        with open(f, "wb") as f:
            _dump_svmlight(X, y, f, one_based, comment, query_id)
