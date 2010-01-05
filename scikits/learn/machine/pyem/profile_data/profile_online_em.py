# /usr/bin/python
# Last Change: Wed Dec 06 08:00 PM 2006 J
import copy

import numpy as N

from gauss_mix import GM
from gmm_em import GMM

def _generate_data(nframes, d, k, mode = 'diag'):
    N.random.seed(0)
    w, mu, va   = GM.gen_param(d, k, mode, spread = 1.5)
    gm          = GM.fromvalues(w, mu, va)
    # Sample nframes frames  from the model
    data        = gm.sample(nframes)

    #++++++++++++++++++++++++++++++++++++++++++
    # Approximate the models with classical EM
    #++++++++++++++++++++++++++++++++++++++++++
    emiter  = 5
    # Init the model
    lgm = GM(d, k, mode)
    gmm = GMM(lgm, 'kmean')
    gmm.init(data)

    gm0    = copy.copy(gmm.gm)
    # The actual EM, with likelihood computation
    like    = N.zeros(emiter)
    for i in range(emiter):
        g, tgd  = gmm.sufficient_statistics(data)
        like[i] = N.sum(N.log(N.sum(tgd, 1)), axis = 0)
        gmm.update_em(data, g)

    return data, gm

nframes = int(5e3)
d       = 1
k       = 2
niter   = 1

def test_v1():
    # Generate test data
    data, gm    = _generate_data(nframes, d, k)
    for i in range(niter):
        iter_1(data, gm)

def test_v2():
    # Generate test data
    data, gm    = _generate_data(nframes, d, k)
    for i in range(niter):
        iter_2(data, gm)

def test_v3():
    # Generate test data
    data, gm    = _generate_data(nframes, d, k)
    for i in range(niter):
        iter_3(data, gm)

def test_v4():
    # Generate test data
    data, gm    = _generate_data(nframes, d, k)
    for i in range(niter):
        iter_4(data, gm)

def iter_1(data, gm):
    """Online EM with original densities + original API"""
    from online_em import OnGMM

    nframes     = data.shape[0]
    ogm         = copy.copy(gm)
    ogmm        = OnGMM(ogm, 'kmean')
    init_data   = data[0:nframes / 20, :]
    ogmm.init(init_data)

    # Forgetting param
    ku		= 0.005
    t0		= 200
    lamb	= 1 - 1/(N.arange(-1, nframes-1) * ku + t0)
    nu0		= 0.2
    nu		= N.zeros((len(lamb), 1))
    nu[0]	= nu0
    for i in range(1, len(lamb)):
        nu[i]	= 1./(1 + lamb[i] / nu[i-1])

    # object version of online EM
    for t in range(nframes):
        ogmm.compute_sufficient_statistics_frame(data[t], nu[t])
        ogmm.update_em_frame()

    ogmm.gm.set_param(ogmm.cw, ogmm.cmu, ogmm.cva)
    print ogmm.cw
    print ogmm.cmu
    print ogmm.cva

def iter_2(data, gm):
    """Online EM with densities2 + original API"""
    from online_em2 import OnGMM

    nframes     = data.shape[0]
    ogm         = copy.copy(gm)
    ogmm        = OnGMM(ogm, 'kmean')
    init_data   = data[0:nframes / 20, :]
    ogmm.init(init_data)

    # Forgetting param
    ku		= 0.005
    t0		= 200
    lamb	= 1 - 1/(N.arange(-1, nframes-1) * ku + t0)
    nu0		= 0.2
    nu		= N.zeros((len(lamb), 1))
    nu[0]	= nu0
    for i in range(1, len(lamb)):
        nu[i]	= 1./(1 + lamb[i] / nu[i-1])

    # object version of online EM
    for t in range(nframes):
        ogmm.compute_sufficient_statistics_frame(data[t], nu[t])
        ogmm.update_em_frame()

    ogmm.gm.set_param(ogmm.cw, ogmm.cmu, ogmm.cva)
    print ogmm.cw
    print ogmm.cmu
    print ogmm.cva

def iter_3(data, gm):
    """Online EM with densities + 1d API"""
    from online_em import OnGMM1d

    #def blop(self, frame, nu):
    #    self.compute_sufficient_statistics_frame(frame, nu)
    #OnGMM.blop  = blop

    nframes     = data.shape[0]
    ogm         = copy.copy(gm)
    ogmm        = OnGMM1d(ogm, 'kmean')
    init_data   = data[0:nframes / 20, :]
    ogmm.init(init_data[:, 0])

    # Forgetting param
    ku		= 0.005
    t0		= 200
    lamb	= 1 - 1/(N.arange(-1, nframes-1) * ku + t0)
    nu0		= 0.2
    nu		= N.zeros((len(lamb), 1))
    nu[0]	= nu0
    for i in range(1, len(lamb)):
        nu[i]	= 1./(1 + lamb[i] / nu[i-1])

    # object version of online EM
    for t in range(nframes):
        #assert ogmm.cw is ogmm.pw
        #assert ogmm.cva is ogmm.pva
        #assert ogmm.cmu is ogmm.pmu
        a, b, c = ogmm.compute_sufficient_statistics_frame(data[t, 0], nu[t])
        ##ogmm.blop(data[t,0], nu[t])
        ogmm.update_em_frame(a, b, c)

    #ogmm.gm.set_param(ogmm.cw, ogmm.cmu, ogmm.cva)
    print ogmm.cw
    print ogmm.cmu
    print ogmm.cva

def iter_4(data, gm):
    """Online EM with densities2 + 1d API"""
    from online_em2 import OnGMM1d

    #def blop(self, frame, nu):
    #    self.compute_sufficient_statistics_frame(frame, nu)
    #OnGMM.blop  = blop

    nframes     = data.shape[0]
    ogm         = copy.copy(gm)
    ogmm        = OnGMM1d(ogm, 'kmean')
    init_data   = data[0:nframes / 20, :]
    ogmm.init(init_data[:, 0])

    # Forgetting param
    ku		= 0.005
    t0		= 200
    lamb	= 1 - 1/(N.arange(-1, nframes-1) * ku + t0)
    nu0		= 0.2
    nu		= N.zeros((len(lamb), 1))
    nu[0]	= nu0
    for i in range(1, len(lamb)):
        nu[i]	= 1./(1 + lamb[i] / nu[i-1])

    # object version of online EM
    def blop():
        #for t in range(nframes):
        #    #assert ogmm.cw is ogmm.pw
        #    #assert ogmm.cva is ogmm.pva
        #    #assert ogmm.cmu is ogmm.pmu
        #    #a, b, c = ogmm.compute_sufficient_statistics_frame(data[t, 0], nu[t])
        #    ###ogmm.blop(data[t,0], nu[t])
        #    #ogmm.update_em_frame(a, b, c)
        #    ogmm.compute_em_frame(data[t, 0], nu[t])
        [ogmm.compute_em_frame(data[t, 0], nu[t]) for t in range(nframes)]
    blop()

    #ogmm.gm.set_param(ogmm.cw, ogmm.cmu, ogmm.cva)
    print ogmm.cw
    print ogmm.cmu
    print ogmm.cva



if __name__ == '__main__':
    #import hotshot, hotshot.stats
    #profile_file    = 'onem1.prof'
    #prof    = hotshot.Profile(profile_file, lineevents=1)
    #prof.runcall(test_v1)
    #p = hotshot.stats.load(profile_file)
    #print p.sort_stats('cumulative').print_stats(20)
    #prof.close()

    #import hotshot, hotshot.stats
    #profile_file    = 'onem2.prof'
    #prof    = hotshot.Profile(profile_file, lineevents=1)
    #prof.runcall(test_v2)
    #p = hotshot.stats.load(profile_file)
    #print p.sort_stats('cumulative').print_stats(20)
    #prof.close()

    import hotshot, hotshot.stats
    profile_file    = 'onem3.prof'
    prof    = hotshot.Profile(profile_file, lineevents=1)
    prof.runcall(test_v3)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()

    import hotshot, hotshot.stats
    profile_file    = 'onem4.prof'
    prof    = hotshot.Profile(profile_file, lineevents=1)
    prof.runcall(test_v4)
    p = hotshot.stats.load(profile_file)
    print p.sort_stats('cumulative').print_stats(20)
    prof.close()
    #test_v1()
    #test_v2()
    #test_v3()
