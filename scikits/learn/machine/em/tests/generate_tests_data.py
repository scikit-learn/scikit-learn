#! /usr/bin/env python
# Last Change: Mon Jun 11 05:00 PM 2007 J

# This script generates some random data used for testing EM implementations.
import copy
import numpy as N
from numpy.testing import set_package_path, restore_path
from scipy.io import savemat, loadmat

set_package_path()
import pyem
restore_path()

from pyem import GM, GMM, EM

def generate_dataset(d, k, mode, nframes):
    """Generate a dataset useful for EM anf GMM testing.
    
    returns:
        data : ndarray
            data from the true model.
        tgm : GM
            the true model (randomly generated)
        gm0 : GM
            the initial model
        gm : GM
            the trained model
    """
    # Generate a model
    w, mu, va = GM.gen_param(d, k, mode, spread = 2.0)
    tgm = GM.fromvalues(w, mu, va)

    # Generate data from the model
    data = tgm.sample(nframes)

    # Run EM on the model, by running the initialization separetely.
    gmm = GMM(GM(d, k, mode), 'test')
    gmm.init_random(data)
    gm0 = copy.copy(gmm.gm)

    gmm = GMM(copy.copy(gmm.gm), 'test')
    em = EM()
    em.train(data, gmm)

    return data, tgm, gm0, gmm.gm

def save_dataset(filename, data, tgm, gm0, gm):
    dic = {'tw': tgm.w, 'tmu': tgm.mu, 'tva': tgm.va,
            'w0': gm0.w, 'mu0' : gm0.mu, 'va0': gm0.va,
            'w': gm.w, 'mu': gm.mu, 'va': gm.va,
            'data': data}
    savemat(filename, dic)

def doall(d, k, mode):
    import pylab as P

    data, tgm, gm0, gm = generate_dataset(d, k, mode, 500)
    filename = mode + '_%dd' % d + '_%dk.mat' % k
    save_dataset(filename, data, tgm, gm0, gm)

    if d == 1:
        P.subplot(2, 1, 1)
        gm0.plot1d()
        h = tgm.plot1d(gpdf = True)
        P.hist(data[:, 0], 20, normed = 1, fill = False)

        P.subplot(2, 1, 2)
        gm.plot1d()
        tgm.plot1d(gpdf = True)
        P.hist(data[:, 0], 20, normed = 1, fill = False)
    else:
        P.subplot(2, 1, 1)
        gm0.plot()
        h = tgm.plot()
        [i.set_color('g') for i in h]
        P.plot(data[:, 0], data[:, 1], '.')

        P.subplot(2, 1, 2)
        gm.plot()
        h = tgm.plot()
        [i.set_color('g') for i in h]
        P.plot(data[:, 0], data[:, 1], '.')

    P.show()

if __name__ == '__main__':
    N.random.seed(0)
    d = 2
    k = 3
    mode = 'full'
    doall(d, k, mode)

    N.random.seed(0)
    d = 2
    k = 3
    mode = 'diag'
    doall(d, k, mode)

    N.random.seed(0)
    d = 1
    k = 4
    mode = 'diag'
    doall(d, k, mode)
