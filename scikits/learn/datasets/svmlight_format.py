
import os.path

import scipy.sparse as sp

from _svmlight_format import _load_svmlight_format

def load_svmlight_format(train_file, test_file=None,
                         n_features=None, buffer_mb=40):

    if not os.path.exists(train_file):
        raise ValueError("Training file doesn't exist")

    ret = []

    data, indices, indptr, labels = _load_svmlight_format(train_file, buffer_mb)

    if n_features is not None:
        shape = (indptr.shape[0] - 1, n_features)
    else:
        shape = None # inferred

    X_train = sp.csr_matrix((data, indices, indptr), shape)

    ret.append(X_train)
    ret.append(labels)

    if test_file is not None:
        if not os.path.exists(test_file):
            raise ValueError("Test file doesn't exist")

        data, indices, indptr, labels = _load_svmlight_format(test_file,
                                                              buffer_mb)

        if n_features is None:
            n_features = X_train.shape[1]

        shape = (indptr.shape[0] - 1, n_features)

        X_test = sp.csr_matrix((data, indices, indptr), shape)

        ret.append(X_test)
        ret.append(labels)

    return tuple(ret)

