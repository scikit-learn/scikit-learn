import numpy as N
from pyem import GM, GMM
import copy

def bench1(mode = 'diag'):
    #===========================================
    # GMM of 20 comp, 20 dimension, 1e4 frames
    #===========================================
    d       = 15
    k       = 30
    nframes = 1e5
    niter   = 10
    mode    = 'diag'

    #+++++++++++++++++++++++++++++++++++++++++++
    # Create an artificial GMM model, samples it
    #+++++++++++++++++++++++++++++++++++++++++++
    print "Generating the mixture"
    # Generate a model with k components, d dimensions
    w, mu, va   = GM.gen_param(d, k, mode, spread = 3)
    # gm          = GM(d, k, mode)
    # gm.set_param(w, mu, va)
    gm          = GM.fromvalues(w, mu, va)

    # Sample nframes frames  from the model
    data    = gm.sample(nframes)

    #++++++++++++++++++++++++
    # Learn the model with EM
    #++++++++++++++++++++++++

    # Init the model
    print "Init a model for learning, with kmean for initialization"
    lgm = GM(d, k, mode)
    gmm = GMM(lgm, 'kmean')
    
    gmm.init(data)
    # Keep the initialized model for drawing
    gm0 = copy.copy(lgm)

    # The actual EM, with likelihood computation
    like    = N.zeros(niter)

    print "computing..."
    for i in range(niter):
        print "iteration %d" % i
        g, tgd  = gmm.sufficient_statistics(data)
        like[i] = N.sum(N.log(N.sum(tgd, 1)))
        gmm.update_em(data, g)

if __name__ == "__main__":
    import profile
    profile.run('bench1()', 'gmmprof')
    import pstats
    p = pstats.Stats('gmmprof')
    print p.sort_stats('cumulative').print_stats(20)

