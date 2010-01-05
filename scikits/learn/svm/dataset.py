__all__ = [
    'LibSvmRegressionDataSet',
    'LibSvmClassificationDataSet',
    'LibSvmOneClassDataSet'
    ]

import numpy as N

import libsvm

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

class LibSvmDataSet:
    def __init__(self, data, kernel):
        self.data = data
        self.kernel = kernel

    def precompute(self):
        n = 1 + len(data) + 1
        a = N.zeros((n,), dtype=svm_node_dtype)
        a[-1] = -1, 0. # end of record marker
        grammat = [a.copy() for i in range(len(data))]
        for i, a in enumerate(grammat):
            a[0] = 0, i + 1 # id
        for i in range(len(data)):
            for j in range(i, len(data)):
                z = self.kernel(data[i], N.transpose(data[j]))
                grammat[i][j+1]['value'] = z
                grammat[j][i+1]['value'] = z
        return LibSvmPrecomputedDataSet(grammat, self.data, self.kernel)

class LibSvmPrecomputedDataSet:
    def __init__(self, grammat, data, kernel):
        self.grammat = grammat
        self.data = data
        self.kernel = kernel

    def extend(self, data):
        raise NotImplementedError

class LibSvmRegressionDataSet(LibSvmDataSet):
    def __init__(self, data):
        f = lambda x: (x[0], convert_to_svm_node(x[1]))
        LibSvmDataSet.__init__(self, map(f, data))

class LibSvmClassificationDataSet(LibSvmDataSet):
    def __init__(self, data):
        labels = N.array(map(lambda x: x[0]), dtype=N.intc)
        assert N.alltrue(labels >= 0), \
            'labels must be non-negative integers'
        f = lambda x: (x[0],convert_to_svm_node(x[1]))
        LibSvmDataSet.__init__(self, map(f, data))

class LibSvmOneClassDataSet(LibSvmDataSet):
    def __init__(self, data):
        f = map(lambda x: tuple([0,convert_to_svm_node(x)]))
        LibSvmDataSet.__init__(self, map(f, data))

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
