# Last Change: Wed Oct 18 06:00 PM 2006 J

import numpy as N
import tables as T

from numpy.random import seed

from gmm_em import multiple_gauss_den
from gauss_mix import GM
from _c_densities import gauss_den

filename    = 'test_mgden.h5';
h5file      = T.openFile(filename, 'w')
h5file.createGroup(h5file.root, 'hyperparams')
h5file.createGroup(h5file.root, 'params')
h5file.createGroup(h5file.root, 'data')

d       = 1
k       = 2
type    = 'diag'
nframes = int(1e3)

h5file.createArray(h5file.root.hyperparams, 'dimension', d)
h5file.createArray(h5file.root.hyperparams, 'type', type)
h5file.createArray(h5file.root.hyperparams, 'nclusters', k)

w, mu, va   = GM.gen_param(d, k, type)

h5file.createArray(h5file.root.params, 'weights', w)
h5file.createArray(h5file.root.params, 'means', mu)
h5file.createArray(h5file.root.params, 'variances', va)

gm      = GM.fromvalues(w, mu, va)
# Sample nframes frames  from the model
data    = gm.sample(nframes)

h5file.createArray(h5file.root.data, 'data', data)

w1, mu1, va1    = GM.gen_param(d, k, type)

out     = multiple_gauss_den(data, mu1, va1)
out1    = gauss_den(data, mu1[0, :], va1[0, :])

h5file.createArray(h5file.root.params, 'w', w1)
h5file.createArray(h5file.root.params, 'mu', mu1)
h5file.createArray(h5file.root.params, 'va', va1)
h5file.createArray(h5file.root.data, 'out', out)

h5file.createArray(h5file.root.params, 'mu1', mu1[0,:])
h5file.createArray(h5file.root.params, 'va1', va1[0,:])
h5file.createArray(h5file.root.data, 'out1', out1)

h5file.close()
