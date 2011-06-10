
import os.path

import scipy.sparse as sp

from _svmlight_format import _load_svmlight_format

def load_svmlight_format(path, buffer_mb=40):
    if not os.path.exists(path):
        raise ValueError("File doesn't exist")
    data, indices, indptr, labels = _load_svmlight_format(path, buffer_mb)
    return sp.csr_matrix((data, indices, indptr)), labels


if __name__ == '__main__':
    import sys
    X, y = load_svmlight_format(sys.argv[1])
    print X.shape
    print y.shape
