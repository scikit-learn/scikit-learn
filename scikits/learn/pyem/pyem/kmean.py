# /usr/bin/python
# Last Change: Mon Aug 28 09:00 PM 2006 J

import numpy as N

def _py_vq(data, code):
    """ Please do not use directly. Use kmean instead"""
    # No attempt to be efficient has been made...
    (n, d)  = data.shape
    (k, d)  = code.shape

    label   = N.zeros(n, int)
    for i in range(n):
        d           = N.sum((data[i, :] - code) ** 2, 1)
        label[i]    = N.argmin(d)

    return label
    
# Try to import pyrex function for vector quantization. If not available,
# falls back on pure python implementation.
#%KMEANIMPORT%
#try:
#    from scipy.cluster.vq import kmeans as kmean
#except ImportError:
#    try:
#        from c_gmm import _vq
#    except:
#        print """c_gmm._vq not found, using pure python implementation instead. 
#        Kmean will be REALLY slow"""
#        _vq = _py_vq
try:
    from sccipy.cluster.vq import vq
    print "using scipy.cluster.vq"
    def _vq(*args, **kw): return vq(*args, **kw)[0]
except ImportError:
    try:
        from c_gmm import _vq
        print "using pyrex vq"
    except ImportError:
        print """c_gmm._vq not found, using pure python implementation instead. 
        Kmean will be REALLY slow"""
        _vq = _py_vq

def kmean(data, init, iter = 10):
    """Simple kmean implementation for EM
    
    returns a tuple (code, label), where code are the final
    centroids, and label are the class label indec for each
    frame (ie row) of data"""

    data    = N.atleast_2d(data)
    init    = N.atleast_2d(init)

    (n, d)  = data.shape
    (k, d1) = init.shape

    if not d == d1:
        msg = "data and init centers do not have same dimensions..."
        raise GmmParamError(msg)
    
    code    = N.asarray(init.copy())
    for i in range(iter):
        # Compute the nearest neighbour for each obs
        # using the current code book
        label   = _vq(data, code)
        # Update the code by computing centroids using the new code book
        for j in range(k):
            code[j,:] = N.mean(data[N.where(label==j)], axis=0) 

    return code, label

# Test functions usable for now
def test_kmean():
    X   = N.array([[3.0, 3], [4, 3], [4, 2],
            [9, 2], [5, 1], [6, 2], [9, 4], 
            [5, 2], [5, 4], [7, 4], [6, 5]])

    initc   = N.concatenate(([[X[0]], [X[1]], [X[2]]])) 

    codet1  = N.array([[3.0000, 3.0000],
            [6.2000, 4.0000], 
            [5.8000, 1.8000]])
            
    codet2  = N.array([[11.0/3, 8.0/3], 
            [6.7500, 4.2500],
            [6.2500, 1.7500]])

    code    = initc.copy()

    code1   = kmean(X, code, 1)[0]
    code2   = kmean(X, code, 2)[0]

    import numpy.testing as testing
    try:
        testing.assert_array_almost_equal(code1, codet1)
        print "kmean test 1 succeded"
    except AssertionError:
        print "kmean test 1 failed"

    try:
        testing.assert_array_almost_equal(code2, codet2)
        print "kmean test 2 succeded"
    except AssertionError:
        print "kmean test 2 failed"

if __name__ == "__main__":
    test_kmean()
