import numpy as np
from scipy.cluster.vq import kmeans2

from gm import GM
from gmm import SStats, Parameters, EM

def initkmeans(data, k):
    # XXX: This is bogus initialization should do better (in kmean with CV)
    (code, label) = kmeans2(data, data[:k], 5, minit='matrix')

    w = np.ones(k) / k
    mu = code.copy()
    va = np.zeros((k, d))
    for c in range(k):
        for i in range(d):
            va[c, i] = np.cov(data[np.where(label==c), i], rowvar=0)

    return w, mu, va

d = 13
k = 15
n = 1e5
print "PB size: %g" % float(n * k * d)

w, mu, va = GM.genparams(d, k)
gm = GM.fromvalues(w, mu, va)

data = gm.sample(n)

# Init the model with kmeans
wi, mui, vai = initkmeans(data, k)

def new_em(data, w, mu, va, niter):
    k = w.size
    k = mu.shape[0]

    params = Parameters.fromvalues(w, mu, va)
    em = EM(hint=1e4)
    em.train(data, params)

    print params.w, params.w.sum()
    return GM.fromvalues(params.w, params.mu, params.va)

def old_em(data, w, mu, va, niter):
    from scikits.learn.machine.em import EM as OEM, GMM as OGMM, GM as OGM

    k = w.size
    k = mu.shape[0]

    lgm = OGM(d, k)
    gmm = OGMM(lgm)

    gmm.gm.w = w.copy()
    gmm.gm.mu = mu.copy()
    gmm.gm.va = va.copy()

    for i in range(niter):
        g, tgd = gmm.compute_responsabilities(data)
        gmm.update_em(data, g)

    return gmm.gm

print "NEW"
newgmm = new_em(data, wi, mui, vai, 10)
print "OLD"
oldgmm = old_em(data, wi, mui, vai, 10)
