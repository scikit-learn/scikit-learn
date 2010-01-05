import numpy as N

from gauss_mix import GM
from gmm_em import GMM, EM

from numpy.random import seed

def test_reg():
    seed(0)
    # Generate data with a few components
    d   = 2
    k   = 1
    n   = 500

    w, mu, va   = GM.gen_param(d, k)
    gm          = GM.fromvalues(w, mu, va)

    data    = gm.sample(n)

    # Try to learn with an insane number of components
    gmm = GMM(GM(d, 30), 'random')

    em  = EM()
    like= em.train(data, gmm, 20, 1e-20)

    # import pylab as P
    # P.subplot(2, 1, 1)
    # P.plot(data[:, 0], data[:, 1], '.')
    # gmm.gm.plot()
    # P.subplot(2, 1, 2)
    # P.plot(like)
    # print like
    # P.show()

if __name__ == "__main__":
    # import hotshot, hotshot.stats
    # profile_file    = 'manyk.prof'
    # prof    = hotshot.Profile(profile_file, lineevents=1)
    # prof.runcall(test_reg)
    # p = hotshot.stats.load(profile_file)
    # print p.sort_stats('cumulative').print_stats(20)
    # prof.close()
    test_reg()

