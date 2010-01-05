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
