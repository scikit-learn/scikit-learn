__all__ = [
    'LibSvmRegressionDataSet',
    'LibSvmClassificationDataSet',
    'LibSvmOneClassDataSet',
    'LibSvmTestDataSet'
    ]

import numpy as N

import libsvm

class LibSvmDataSet:
    def __init__(self, data):
        self.data = data

    def getgamma(self):
        maxlen = 0
        for y, x in self.data:
            maxlen = N.maximum(maxlen, x['index'].max())
        return 1.0 / maxlen
    gamma = property(getgamma, 'Gamma parameter for RBF kernel')

    def precompute(self):
        return LibSvmPrecomputedDataSet(self.data)

class LibSvmPrecomputedDataSet:
    def __init__(self, origdata):
        # XXX look at using a list of vectors instead of a matrix when
        # the size of the precomputed dataset gets huge. This should
        # avoid problems with heap fragmentation, especially on
        # Windows.

        n = len(origdata)
        # extra columns for id and end of record marker
        grammat = N.empty((n, n+2), dtype=libsvm.svm_node_dtype)
        # calculate Gram matrix
        for i, (y1, x1) in enumerate(origdata):
            # set id and end of record fields
            grammat[i,0], grammat[i,-1] = (0, i), (-1, 0.0)
            for j, (y2, x2) in enumerate(origdata[i:]):
                # Gram matrix is symmetric, so calculate dot product
                # once and store it in both required locations
                z = svm_node_dot(x1, x2)
                # fix index so we assign to the right place
                j += i
                grammat[i, j+1]['value'] = z
                grammat[j, i+1]['value'] = z
        self.grammat = grammat
        self.data = zip(map(lambda x: x[0], origdata), grammat)

class LibSvmRegressionDataSet(LibSvmDataSet):
    def __init__(self, origdata):
        data = map(lambda x: (x[0], convert_to_svm_node(x[1])), origdata)
        LibSvmDataSet.__init__(self, data)

class LibSvmClassificationDataSet(LibSvmDataSet):
    def __init__(self, origdata):
        labels = N.array(map(lambda x: x[0], origdata), dtype=N.intc)
        labels.sort()
        self.labels = labels

        data = map(lambda x: (x[0],convert_to_svm_node(x[1])), origdata)
        LibSvmDataSet.__init__(self, data)

class LibSvmOneClassDataSet(LibSvmDataSet):
    def __init__(self, origdata):
        data = map(lambda x: tuple([0,convert_to_svm_node(x)]), origdata)
        LibSvmDataSet.__init__(self, data)

class LibSvmTestDataSet:
    def __init__(self, origdata):
        self.data = map(lambda x: convert_to_svm_node(x), origdata)

def convert_to_svm_node(x):
    y = N.empty(len(x)+1, dtype=libsvm.svm_node_dtype)
    y[-1] = (-1, 0.)
    if isinstance(x, dict):
        x = x.items()
    if isinstance(x, list):
        x.sort(cmp=lambda x,y: cmp(x[0],y[0]))
        y[:-1] = x
    else:
        y['index'][:-1] = N.arange(1,len(x)+1)
        y['value'][:-1] = x
    assert N.alltrue(y[:-1]['index'] >= 1), \
        'indexes must be positive'
    assert len(x) == len(N.unique(y[:-1]['index'])), \
        'indexes must be unique'
    return y

def svm_node_dot(x, y):
    # associate node indexes with array indexes
    xidx = dict(zip(x['index'][:-1],range(0,len(x))))
    yidx = dict(zip(y['index'][:-1],range(0,len(y))))
    # indexes in either vector
    indexes = N.unique(N.hstack([x['index'],y['index']]))
    z = 0.
    for j in indexes:
        if j in xidx and j in yidx:
            # dot if index is present in both vectors
            z += x['value'][xidx[j]]*y['value'][yidx[j]]
    return z
