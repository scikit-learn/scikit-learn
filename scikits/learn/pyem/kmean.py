# /usr/bin/python
# Last Change: Thu Sep 28 01:00 PM 2006 J

#TODO:
#   - a demo for kmeans

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
    from scipy.cluster.vq import vq
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
    """Simple kmean implementation for EM. Runs iter iterations.
    
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

if __name__ == "__main__":
    pass
