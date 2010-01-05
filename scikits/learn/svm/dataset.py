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
        self.iddatamap = {}
        for y, x in data:
            key = x.__array_interface__['data'][0]
            self.iddatamap[key] = x
        self.__len__ = self.data.__len__
        self.__iter__ = self.data.__iter__

    def getgamma(self):
        maxlen = 0
        for y, x in self.data:
            maxlen = N.maximum(maxlen, x['index'].max())
        return 1.0 / maxlen
    gamma = property(getgamma, 'Gamma parameter for RBF kernel')

    def precompute(self, kernel):
        return LibSvmPrecomputedDataSet(kernel, self.data)

    def _create_svm_problem(self):
        return libsvm.create_svm_problem(self.data)

    def _update_svm_parameter(self, param):
        # XXX we can handle gamma=None here
        pass

class LibSvmPrecomputedDataSet:
    def __init__(self, kernel, origdata=None):
        self.kernel = kernel
        self.origdata = origdata
        if origdata is None: return

        self.iddatamap = {}

        # Create Gram matrix as a list of vectors which an extra entry
        # for the id field.
        n = len(origdata)
        grammat = [N.empty((n + 1,), dtype=libsvm.svm_node_dtype)
                   for i in range(n)]
        self.grammat = grammat

        # Calculate Gram matrix. Refer to Kernel::kernel_precomputed
        # in svm.cpp to see how this precomputed setup works.
        for i, (yi, xi) in enumerate(origdata):
            id = i + 1
            grammat[i][0] = 0, id
            # Map id to original vector so that we can find it again
            # after the model has been trained. libsvm essentially
            # provides the ids of the support vectors.
            self.iddatamap[id] = xi
            for j, (yj, xj) in enumerate(origdata[i:]):
                # Gram matrix is symmetric, so calculate dot product
                # once and store it in both required locations
                z = self.kernel(xi, xj, svm_node_dot)
                # fix index so we assign to the right place
                j += i
                grammat[i][j + 1] = 0, z
                grammat[j][i + 1] = 0, z

    def __len__(self):
        return len(self.origdata)

    def __getitem__(self, key):
        return self.iddatamap[key]

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
        n = len(self.origdata) + len(dataset.data) + 1
        newgrammat = []

        # copy original Gram matrix
        for i in range(len(self.origdata)):
            newrow = N.zeros((n,), dtype=libsvm.svm_node_dtype)
            oldrow = self.grammat[i]
            newrow[:len(oldrow)] = oldrow
            newgrammat.append(newrow)

        # prepare Gram matrix for new data
        for i in range(len(dataset.data)):
            row = N.zeros((n,), dtype=libsvm.svm_node_dtype)
            newgrammat.append(row)

        newiddatamap = dict(self.iddatamap.items())
        m = len(self.origdata)
        for i, (yi, xi) in enumerate(dataset.data):
            i += m
            for j, (yj, xj) in enumerate(self.origdata):
                z = self.kernel(xi, xj, svm_node_dot)
                newgrammat[i][j + 1] = 0, z
                newgrammat[j][i + 1] = 0, z
        for i, (yi, xi) in enumerate(dataset.data):
            k = m + i
            id = k + 1
            newgrammat[k][0] = 0, id
            newiddatamap[id] = xi
            for j, (yj, xj) in enumerate(dataset.data[i:]):
                z = self.kernel(xi, xj, svm_node_dot)
                j += k
                newgrammat[k][j + 1] = 0, z
                newgrammat[j][k + 1] = 0, z

        newdataset = self.__class__(self.kernel)
        newdataset.origdata = self.origdata + dataset.data
        newdataset.iddatamap = newiddatamap
        newdataset.grammat = newgrammat
        return newdataset

    def _create_svm_problem(self):
        return libsvm.create_svm_problem(self.data)

    def _update_svm_parameter(self, param):
        param.kernel_type = libsvm.PRECOMPUTED

class LibSvmRegressionDataSet(LibSvmDataSet):
    def __init__(self, y, x):
        origdata = zip(y, x)
        data = [(x[0], convert_to_svm_node(x[1])) for x in origdata]
        LibSvmDataSet.__init__(self, data)

class LibSvmClassificationDataSet(LibSvmDataSet):
    def __init__(self, labels, x):
        origdata = zip(labels, x)
        data = [(x[0], convert_to_svm_node(x[1])) for x in origdata]
        LibSvmDataSet.__init__(self, data)

class LibSvmOneClassDataSet(LibSvmDataSet):
    def __init__(self, x):
        data = [(0, convert_to_svm_node(y)) for y in x]
        LibSvmDataSet.__init__(self, data)

class LibSvmTestDataSet:
    def __init__(self, origdata):
        self.data = map(lambda x: convert_to_svm_node(x), origdata)
        self.__len__ = self.data.__len__
        self.__iter__ = self.data.__iter__
        self.__getitem__ = self.data.__getitem__

def convert_to_svm_node(x):
    y = N.empty(len(x) + 1, dtype=libsvm.svm_node_dtype)
    y[-1] = -1, 0.
    if isinstance(x, dict):
        x = x.items()
    if isinstance(x, list):
        x.sort(cmp=lambda x,y: cmp(x[0],y[0]))
        y[:-1] = x
    else:
        y['index'][:-1] = N.arange(1,len(x) + 1)
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
            z += x['value'][xidx[j]] * y['value'][yidx[j]]
    return z
