from ctypes import c_double, POINTER, cast
import numpy as N

import libsvm

__all__ = [
    'LibSvmRegressionDataSet',
    'LibSvmClassificationDataSet',
    'LibSvmOneClassDataSet',
    'LibSvmTestDataSet'
    ]

class LibSvmDataSet:
    def __init__(self, data):
        self.data = data

    def getgamma(self):
        maxlen = 0
        for y, x in self.data:
            maxlen = N.maximum(maxlen, x['index'].max())
        return 1.0 / maxlen
    gamma = property(getgamma, 'Gamma parameter for RBF kernel')

    def precompute(self, kernel):
        return LibSvmPrecomputedDataSet(kernel, self.data)

    def create_svm_problem(self):
        problem = libsvm.svm_problem()
        problem.l = len(self.data)
        y = (c_double*problem.l)()
        x = (POINTER(libsvm.svm_node)*problem.l)()
        for i, (yi, xi) in enumerate(self.data):
            y[i] = yi
            x[i] = cast(xi.ctypes.data, POINTER(libsvm.svm_node))
        problem.x = x
        problem.y = y
        return problem

class LibSvmPrecomputedDataSet:
    def __init__(self, kernel, origdata=None):
        self.kernel = kernel
        self.origdata = origdata
        if origdata is None: return

        self.iddatamap = {}

        # Create Gram matrix as a list of vectors that have extra
        # entries for id and end of record marker.
        n = len(origdata)
        grammat = [N.empty((n+2,), dtype=libsvm.svm_node_dtype)
                   for i in range(n)]
        self.grammat = grammat

        # Calculate Gram matrix. Refer to Kernel::kernel_precomputed
        # in svm.cpp to see how this precomputed setup works.
        for i, (y1, x1) in enumerate(origdata):
            id = i + 1
            # XXX possible numpy bug
            #grammat[i][[0,-1]] = (0, id), (-1, 0.0)
            grammat[i][0] = 0, id
            grammat[i][-1] = -1, 0.0
            for j, (y2, x2) in enumerate(origdata[i:]):
                # Gram matrix is symmetric, so calculate dot product
                # once and store it in both required locations
                z = kernel(x1, x2, svm_node_dot)
                # fix index so we assign to the right place
                j += i
                grammat[i][j+1] = 0, z
                grammat[j][i+1] = 0, z
            # Map id to original vector so that we can find it again
            # after the model has been trained. libsvm essentially
            # provides the ids of the support vectors.
            self.iddatamap[id] = x1
    
    def getdata(self):
        return zip(map(lambda x: x[0], self.origdata), self.grammat)
    data = property(getdata)

    def combine_inplace(self, dataset):
        """
        Combine this dataset with another dataset by calculating the
        new part of the Gram matrix in place.
        """
        # XXX N.resize is our friend here
        raise NotImplementedError

    def combine(self, dataset):
        """
        Combine this dataset with another dataset by extending the
        Gram matrix with the new inner products into a new matrix.
        """
        n = len(self.origdata) + len(dataset.data)
        newgrammat = []

        # copy original Gram matrix
        for i in range(len(self.origdata)):
            row = N.empty((n,), dtype=libsvm.svm_node_dtype)
            row[:-1] = self.grammat[i]
            newgrammat.append(row)

        # copy id->vector map
        newiddatamap = dict(self.iddatamap.items())

        # prepare Gram matrix for new data
        for i in range(len(dataset.data)):
            id = i + len(self.origdata) + 1
            row = N.empty((n,), dtype=libsvm.svm_node_dtype)
            row[[0,-1]] = (0, id), (-1, 0.0)
            newgrammat.append(row)
            newiddatamap[id] = dataset.data[i][1]

        newdataset = self.__class__(self.kernel)
        newdataset.origdata = self.origdata + dataset.data
        newdataset.iddatamap = newiddatamap
        newdataset.grammat = newgrammat

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
    y[-1] = -1, 0.
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
