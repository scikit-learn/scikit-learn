import numpy as np

from numpy.testing import assert_array_equal, assert_array_almost_equal
from unittest import TestCase

from sklearn.datasets.samples_generator import make_spd_matrix
from sklearn import hmm
from sklearn.utils.extmath import logsumexp


np.seterr(all='warn')


class SeedRandomNumberGeneratorTestCase(TestCase):
    seed = 9

    def __init__(self, *args, **kwargs):
        self.setUp()
        TestCase.__init__(self, *args, **kwargs)

    def setUp(self):
        self.prng = np.random.RandomState(self.seed)


class TestBaseHMM(SeedRandomNumberGeneratorTestCase):

    class StubHMM(hmm._BaseHMM):

        def _compute_log_likelihood(self, X):
            return self.framelogprob

        def _generate_sample_from_state(self):
            pass

        def _init(self):
            pass

    def test_prune_states_no_pruning(self):
        h = self.StubHMM(10)
        lattice_frame = np.arange(h.n_components)

        idx = h._prune_states(lattice_frame, None, - np.Inf)
        assert_array_equal(idx, range(h.n_components))

    def test_prune_states_rank(self):
        h = self.StubHMM(10)
        lattice_frame = np.arange(h.n_components)

        idx = h._prune_states(lattice_frame, 1, -np.Inf)
        assert_array_equal(idx, [lattice_frame.argmax()])

    def test_prune_states_beam(self):
        h = self.StubHMM(10)
        lattice_frame = np.arange(h.n_components)

        beamlogprob = -h.n_components / 2
        idx = h._prune_states(lattice_frame, None, beamlogprob)
        refidx, = np.nonzero(lattice_frame >= -beamlogprob)
        assert_array_equal(idx, refidx)

    def setup_example_hmm(self):
        # Example from http://en.wikipedia.org/wiki/Forward-backward_algorithm
        h = self.StubHMM(2)
        h.transmat = [[0.7, 0.3], [0.3, 0.7]]
        h.start_prob = [0.5, 0.5]
        framelogprob = np.log([[0.9, 0.2],
                               [0.9, 0.2],
                               [0.1, 0.8],
                               [0.9, 0.2],
                               [0.9, 0.2]])
        # Add dummy observations to stub.
        h.framelogprob = framelogprob
        return h, framelogprob

    def test_init(self):
        h, framelogprob = self.setup_example_hmm()

        for params in [('transmat',), ('startprob', 'transmat')]:
            d = dict((x, getattr(h, x)) for x in params)
            h2 = self.StubHMM(h.n_components, **d)
            self.assertEqual(h.n_components, h2.n_components)
            for p in params:
                assert_array_almost_equal(getattr(h, p), getattr(h2, p))

    def test_do_forward_pass(self):
        h, framelogprob = self.setup_example_hmm()

        logprob, fwdlattice = h._do_forward_pass(framelogprob)

        reflogprob = -3.3725
        self.assertAlmostEqual(logprob, reflogprob, places=4)

        reffwdlattice = np.array([[0.4500, 0.1000],
                                  [0.3105, 0.0410],
                                  [0.0230, 0.0975],
                                  [0.0408, 0.0150],
                                  [0.0298, 0.0046]])
        assert_array_almost_equal(np.exp(fwdlattice), reffwdlattice, 4)

    def test_do_backward_pass(self):
        h, framelogprob = self.setup_example_hmm()

        fakefwdlattice = np.zeros((len(framelogprob), 2))
        bwdlattice = h._do_backward_pass(framelogprob, fakefwdlattice)

        refbwdlattice = np.array([[0.0661, 0.0455],
                                  [0.0906, 0.1503],
                                  [0.4593, 0.2437],
                                  [0.6900, 0.4100],
                                  [1.0000, 1.0000]])
        assert_array_almost_equal(np.exp(bwdlattice), refbwdlattice, 4)

    def test_do_viterbi_pass(self):
        h, framelogprob = self.setup_example_hmm()

        logprob, state_sequence = h._do_viterbi_pass(framelogprob)

        refstate_sequence = [0, 0, 1, 0, 0]
        assert_array_equal(state_sequence, refstate_sequence)

        reflogprob = -4.4590
        self.assertAlmostEqual(logprob, reflogprob, places=4)

    def test_eval(self):
        h, framelogprob = self.setup_example_hmm()
        nobs = len(framelogprob)

        logprob, posteriors = h.eval([])

        assert_array_almost_equal(posteriors.sum(axis=1), np.ones(nobs))

        reflogprob = -3.3725
        self.assertAlmostEqual(logprob, reflogprob, places=4)

        refposteriors = np.array([[0.8673, 0.1327],
                                  [0.8204, 0.1796],
                                  [0.3075, 0.6925],
                                  [0.8204, 0.1796],
                                  [0.8673, 0.1327]])
        assert_array_almost_equal(posteriors, refposteriors, decimal=4)

    def test_hmm_eval_consistent_with_gmm(self):
        n_components = 8
        nobs = 10
        h = self.StubHMM(n_components)

        # Add dummy observations to stub.
        framelogprob = np.log(self.prng.rand(nobs, n_components))
        h.framelogprob = framelogprob

        # If startprob and transmat are uniform across all states (the
        # default), the transitions are uninformative - the model
        # reduces to a GMM with uniform mixing weights (in terms of
        # posteriors, not likelihoods).
        logprob, hmmposteriors = h.eval([])

        assert_array_almost_equal(hmmposteriors.sum(axis=1), np.ones(nobs))

        norm = logsumexp(framelogprob, axis=1)[:, np.newaxis]
        gmmposteriors = np.exp(framelogprob - np.tile(norm, (1, n_components)))
        assert_array_almost_equal(hmmposteriors, gmmposteriors)

    def test_hmm_decode_consistent_with_gmm(self):
        n_components = 8
        nobs = 10
        h = self.StubHMM(n_components)

        # Add dummy observations to stub.
        framelogprob = np.log(self.prng.rand(nobs, n_components))
        h.framelogprob = framelogprob

        # If startprob and transmat are uniform across all states (the
        # default), the transitions are uninformative - the model
        # reduces to a GMM with uniform mixing weights (in terms of
        # posteriors, not likelihoods).
        viterbi_ll, state_sequence = h.decode([])

        norm = logsumexp(framelogprob, axis=1)[:, np.newaxis]
        gmmposteriors = np.exp(framelogprob - np.tile(norm, (1, n_components)))
        gmmstate_sequence = gmmposteriors.argmax(axis=1)
        assert_array_equal(state_sequence, gmmstate_sequence)

    def test_base_hmm_attributes(self):
        n_components = 20
        startprob = self.prng.rand(n_components)
        startprob = startprob / startprob.sum()
        transmat = prng.rand(n_components, n_components)
        transmat /= np.tile(transmat.sum(axis=1)
                [:, np.newaxis], (1, n_components))

        h = self.StubHMM(n_components)

        self.assertEquals(h.n_components, n_components)

        h.startprob = startprob
        assert_array_almost_equal(h.startprob, startprob)

        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          2 * startprob)
        self.assertRaises(ValueError, h.__setattr__, 'startprob', [])
        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          np.zeros((n_components - 2, 2)))

        h.transmat = transmat
        assert_array_almost_equal(h.transmat, transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          2 * transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat', [])
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          np.zeros((n_components - 2, n_components)))


def train_hmm_and_keep_track_of_log_likelihood(hmm, obs, n_iter=1, **kwargs):
    hmm.fit(obs, n_iter=1, **kwargs)
    loglikelihoods = []
    for n in xrange(n_iter):
        hmm.fit(obs, n_iter=1, init_params='', **kwargs)
        loglikelihoods.append(sum(hmm.score(x) for x in obs))
    return loglikelihoods

prng = np.random.RandomState(10)


class GaussianHMMParams(object):
    n_components = 3
    n_features = 3
    startprob = prng.rand(n_components)
    startprob = startprob / startprob.sum()
    transmat = np.random.rand(n_components, n_components)
    transmat /= np.tile(transmat.sum(axis=1)[:, np.newaxis], (1, n_components))
    means = prng.randint(-20, 20, (n_components, n_features))
    covars = {'spherical': (1.0 + 2 * prng.rand(n_components)) ** 2,
              'tied': (make_spd_matrix(n_features, random_state=0)
                  + np.eye(n_features)),
              'diag': (1.0 + 2 * prng.rand(n_components, n_features)) ** 2,
              'full': np.array(
                  [make_spd_matrix(n_features,
                      random_state=0) + np.eye(n_features)
                   for x in xrange(n_components)])}
    expanded_covars = {'spherical': [np.eye(n_features) * cov
                                     for cov in covars['spherical']],
                       'diag': [np.diag(cov) for cov in covars['diag']],
                       'tied': [covars['tied']] * n_components,
                       'full': covars['full']}


class GaussianHMMTester(GaussianHMMParams):

    def test_bad_cvtype(self):
        hmm.GaussianHMM(20, self.cvtype)
        self.assertRaises(ValueError, hmm.GaussianHMM, 20, 'badcvtype')

    def test_attributes(self):
        h = hmm.GaussianHMM(self.n_components, self.cvtype)

        self.assertEquals(h.n_components, self.n_components)
        self.assertEquals(h.cvtype, self.cvtype)

        h.startprob = self.startprob
        assert_array_almost_equal(h.startprob, self.startprob)
        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          2 * self.startprob)
        self.assertRaises(ValueError, h.__setattr__, 'startprob', [])
        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          np.zeros((self.n_components - 2, self.n_features)))

        h.transmat = self.transmat
        assert_array_almost_equal(h.transmat, self.transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          2 * self.transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat', [])
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          np.zeros((self.n_components - 2, self.n_components)))

        h.means = self.means
        assert_array_almost_equal(h.means, self.means)
        self.assertEquals(h.n_features, self.n_features)
        self.assertRaises(ValueError, h.__setattr__, 'means', [])
        self.assertRaises(ValueError, h.__setattr__, 'means',
                          np.zeros((self.n_components - 2, self.n_features)))

        h.covars = self.covars[self.cvtype]
        assert_array_almost_equal(h.covars, self.expanded_covars[self.cvtype])
        #self.assertRaises(ValueError, h.__setattr__, 'covars', [])
        #self.assertRaises(ValueError, h.__setattr__, 'covars',
        #                  np.zeros((self.n_components - 2, self.n_features)))

    def test_eval_and_decode(self):
        h = hmm.GaussianHMM(self.n_components, self.cvtype)
        h.means = self.means
        h.covars = self.covars[self.cvtype]

        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        h.means = 20 * h.means

        gaussidx = np.repeat(range(self.n_components), 5)
        nobs = len(gaussidx)
        obs = self.prng.randn(nobs, self.n_features) + h.means[gaussidx]

        ll, posteriors = h.eval(obs)

        self.assertEqual(posteriors.shape, (nobs, self.n_components))
        assert_array_almost_equal(posteriors.sum(axis=1), np.ones(nobs))

        viterbi_ll, stateseq = h.decode(obs)
        assert_array_equal(stateseq, gaussidx)

    def test_rvs(self, n=1000):
        h = hmm.GaussianHMM(self.n_components, self.cvtype)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        h.means = 20 * self.means
        h.covars = np.maximum(self.covars[self.cvtype], 0.1)
        h.startprob = self.startprob

        samples = h.rvs(n)
        self.assertEquals(samples.shape, (n, self.n_features))

    def test_fit(self, params='stmc', n_iter=25, verbose=False, **kwargs):
        np.random.seed(0)
        h = hmm.GaussianHMM(self.n_components, self.cvtype)
        h.startprob = self.startprob
        h.transmat = hmm.normalize(self.transmat
                + np.diag(self.prng.rand(self.n_components)), 1)
        h.means = 20 * self.means
        h.covars = self.covars[self.cvtype]

        # Create training data by sampling from the HMM.
        train_obs = [h.rvs(n=10) for x in xrange(10)]

        # Mess up the parameters and see if we can re-learn them.
        h.fit(train_obs, n_iter=0)

        trainll = train_hmm_and_keep_track_of_log_likelihood(
            h, train_obs, n_iter=n_iter, params=params, **kwargs)[1:]
        if not np.all(np.diff(trainll) > 0) and verbose:
            print
            print ('Test train: %s (%s)\n  %s\n  %s'
                   % (self.cvtype, params, trainll, np.diff(trainll)))
        delta_min = np.diff(trainll).min()
        self.assertTrue(
            delta_min > -0.8,
            "The min nll increase is %f which is lower than the admissible"
            " threshold of %f, for model %s. The likelihoods are %s."
                % (delta_min, -0.8, self.cvtype, trainll))

    def test_fit_works_on_sequences_of_different_length(self):
        obs = [self.prng.rand(3, self.n_features),
               self.prng.rand(4, self.n_features),
               self.prng.rand(5, self.n_features)]

        h = hmm.GaussianHMM(self.n_components, self.cvtype)
        # This shouldn't raise
        # ValueError: setting an array element with a sequence.
        h.fit(obs)

    def test_fit_with_priors(self, params='stmc', n_iter=10,
                             verbose=False):
        startprob_prior = 10 * self.startprob + 2.0
        transmat_prior = 10 * self.transmat + 2.0
        means_prior = self.means
        means_weight = 2.0
        covars_weight = 2.0
        if self.cvtype in ('full', 'tied'):
            covars_weight += self.n_features
        covars_prior = self.covars[self.cvtype]

        h = hmm.GaussianHMM(self.n_components, self.cvtype)
        h.startprob = self.startprob
        h.startprob_prior = startprob_prior
        h.transmat = hmm.normalize(self.transmat
                + np.diag(self.prng.rand(self.n_components)), 1)
        h.transmat_prior = transmat_prior
        h.means = 20 * self.means
        h.means_prior = means_prior
        h.means_weight = means_weight
        h.covars = self.covars[self.cvtype]
        h.covars_prior = covars_prior
        h.covars_weight = covars_weight

        # Create training data by sampling from the HMM.
        train_obs = [h.rvs(n=10) for x in xrange(10)]

        # Mess up the parameters and see if we can re-learn them.
        h.fit(train_obs[:1], n_iter=0)

        trainll = train_hmm_and_keep_track_of_log_likelihood(
            h, train_obs, n_iter=n_iter, params=params)[1:]
        if not np.all(np.diff(trainll) > 0) and verbose:
            print
            print ('Test MAP train: %s (%s)\n  %s\n  %s'
                   % (self.cvtype, params, trainll, np.diff(trainll)))
        self.assertTrue(np.all(np.diff(trainll) > -0.5))


class TestGaussianHMMWithSphericalCovars(GaussianHMMTester,
                                         SeedRandomNumberGeneratorTestCase):
    cvtype = 'spherical'

    def test_fit_startprob_and_transmat(self):
        self.test_fit('st')


class TestGaussianHMMWithDiagonalCovars(GaussianHMMTester,
                                        SeedRandomNumberGeneratorTestCase):
    cvtype = 'diag'


class TestGaussianHMMWithTiedCovars(GaussianHMMTester,
                                    SeedRandomNumberGeneratorTestCase):
    cvtype = 'tied'


class TestGaussianHMMWithFullCovars(GaussianHMMTester,
                                    SeedRandomNumberGeneratorTestCase):
    cvtype = 'full'


class MultinomialHMMParams(object):
    """Using example from http://en.wikipedia.org/wiki/Hidden_Markov_model
    and http://en.wikipedia.org/wiki/Viterbi_algorithm"""
    n_components = 2   # ('Rainy', 'Sunny')
    n_symbols = 3  # ('walk', 'shop', 'clean')
    emissionprob = [[0.1, 0.4, 0.5], [0.6, 0.3, 0.1]]
    startprob = [0.6, 0.4]
    transmat = [[0.7, 0.3], [0.4, 0.6]]


class TestMultinomialHMM(MultinomialHMMParams,
                         SeedRandomNumberGeneratorTestCase):

    def test_wikipedia_viterbi_example(self):
        # From http://en.wikipedia.org/wiki/Viterbi_algorithm:
        # "This reveals that the observations ['walk', 'shop', 'clean']
        # were most likely generated by states ['Sunny', 'Rainy',
        # 'Rainy'], with probability 0.01344."
        observations = [0, 1, 2]

        h = hmm.MultinomialHMM(self.n_components,
                               startprob=self.startprob,
                               transmat=self.transmat)
        h.emissionprob = self.emissionprob
        logprob, state_sequence = h.decode(observations)

        self.assertAlmostEqual(np.exp(logprob), 0.01344)
        assert_array_equal(state_sequence, [1, 0, 0])

    def test_attributes(self):
        h = hmm.MultinomialHMM(self.n_components)

        self.assertEquals(h.n_components, self.n_components)

        h.startprob = self.startprob
        assert_array_almost_equal(h.startprob, self.startprob)
        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          2 * self.startprob)
        self.assertRaises(ValueError, h.__setattr__, 'startprob', [])
        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          np.zeros((self.n_components - 2, self.n_symbols)))

        h.transmat = self.transmat
        assert_array_almost_equal(h.transmat, self.transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          2 * self.transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat', [])
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          np.zeros((self.n_components - 2, self.n_components)))

        h.emissionprob = self.emissionprob
        assert_array_almost_equal(h.emissionprob, self.emissionprob)
        self.assertRaises(ValueError, h.__setattr__, 'emissionprob', [])
        self.assertRaises(ValueError, h.__setattr__, 'emissionprob',
                          np.zeros((self.n_components - 2, self.n_symbols)))
        self.assertEquals(h.n_symbols, self.n_symbols)

    def test_eval(self):
        h = hmm.MultinomialHMM(self.n_components,
                               startprob=self.startprob,
                               transmat=self.transmat)
        h.emissionprob = self.emissionprob
        idx = np.repeat(range(self.n_components), 10)
        nobs = len(idx)
        obs = [int(x) for x in np.floor(self.prng.rand(nobs) * self.n_symbols)]

        ll, posteriors = h.eval(obs)

        self.assertEqual(posteriors.shape, (nobs, self.n_components))
        assert_array_almost_equal(posteriors.sum(axis=1), np.ones(nobs))

    def test_rvs(self, n=1000):
        h = hmm.MultinomialHMM(self.n_components,
                               startprob=self.startprob,
                               transmat=self.transmat)
        h.emissionprob = self.emissionprob
        samples = h.rvs(n)
        self.assertEquals(len(samples), n)
        self.assertEquals(len(np.unique(samples)), self.n_symbols)

    def test_fit(self, params='ste', n_iter=15, verbose=False, **kwargs):
        np.random.seed(0)
        h = hmm.MultinomialHMM(self.n_components,
                               startprob=self.startprob,
                               transmat=self.transmat)
        h.emissionprob = self.emissionprob

        # Create training data by sampling from the HMM.
        train_obs = [h.rvs(n=10) for x in xrange(10)]

        # Mess up the parameters and see if we can re-learn them.
        h.startprob = hmm.normalize(self.prng.rand(self.n_components))
        h.transmat = hmm.normalize(self.prng.rand(self.n_components,
                                                  self.n_components), axis=1)
        h.emissionprob = hmm.normalize(
            self.prng.rand(self.n_components, self.n_symbols), axis=1)

        trainll = train_hmm_and_keep_track_of_log_likelihood(
            h, train_obs, n_iter=n_iter, params=params, **kwargs)[1:]
        if not np.all(np.diff(trainll) > 0) and verbose:
            print
            print 'Test train: (%s)\n  %s\n  %s' % (params, trainll,
                                                    np.diff(trainll))
        self.assertTrue(np.all(np.diff(trainll) > 0))

    def test_fit_emissionprob(self):
        self.test_fit('e')


def create_random_gmm(n_mix, n_features, cvtype, prng=prng):
    from sklearn import mixture

    g = mixture.GMM(n_mix, cvtype=cvtype)
    g.means = prng.randint(-20, 20, (n_mix, n_features))
    mincv = 0.1
    g.covars = {
        'spherical': (mincv + mincv * prng.rand(n_mix)) ** 2,
        'tied': (make_spd_matrix(n_features, random_state=prng)
                 + mincv * np.eye(n_features)),
        'diag': (mincv + mincv * prng.rand(n_mix, n_features)) ** 2,
        'full': np.array(
            [make_spd_matrix(n_features, random_state=prng)
             + mincv * np.eye(n_features) for x in xrange(n_mix)])
    }[cvtype]
    g.weights = hmm.normalize(prng.rand(n_mix))
    return g


class GMMHMMParams(object):
    n_components = 3
    n_mix = 2
    n_features = 2
    cvtype = 'diag'
    startprob = prng.rand(n_components)
    startprob = startprob / startprob.sum()
    transmat = prng.rand(n_components, n_components)
    transmat /= np.tile(transmat.sum(axis=1)[:, np.newaxis], (1, n_components))


class TestGMMHMM(GMMHMMParams, SeedRandomNumberGeneratorTestCase):

    def setUp(self):
        self.prng = np.random.RandomState(self.seed)
        self.gmms = []
        for state in xrange(self.n_components):
            self.gmms.append(create_random_gmm(
                self.n_mix, self.n_features, self.cvtype, prng=self.prng))

    def test_attributes(self):
        h = hmm.GMMHMM(self.n_components, cvtype=self.cvtype)

        self.assertEquals(h.n_components, self.n_components)

        h.startprob = self.startprob
        assert_array_almost_equal(h.startprob, self.startprob)
        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          2 * self.startprob)
        self.assertRaises(ValueError, h.__setattr__, 'startprob', [])
        self.assertRaises(ValueError, h.__setattr__, 'startprob',
                          np.zeros((self.n_components - 2, self.n_features)))

        h.transmat = self.transmat
        assert_array_almost_equal(h.transmat, self.transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          2 * self.transmat)
        self.assertRaises(ValueError, h.__setattr__, 'transmat', [])
        self.assertRaises(ValueError, h.__setattr__, 'transmat',
                          np.zeros((self.n_components - 2, self.n_components)))

    def test_eval_and_decode(self):
        h = hmm.GMMHMM(self.n_components, gmms=self.gmms)
        # Make sure the means are far apart so posteriors.argmax()
        # picks the actual component used to generate the observations.
        for g in h.gmms:
            g.means *= 20

        refstateseq = np.repeat(range(self.n_components), 5)
        nobs = len(refstateseq)
        obs = [h.gmms[x].rvs(1).flatten() for x in refstateseq]

        ll, posteriors = h.eval(obs)

        self.assertEqual(posteriors.shape, (nobs, self.n_components))
        assert_array_almost_equal(posteriors.sum(axis=1), np.ones(nobs))

        viterbi_ll, stateseq = h.decode(obs)
        assert_array_equal(stateseq, refstateseq)

    def test_rvs(self, n=1000):
        h = hmm.GMMHMM(self.n_components, self.cvtype,
                       startprob=self.startprob, transmat=self.transmat,
                       gmms=self.gmms)
        samples = h.rvs(n)
        self.assertEquals(samples.shape, (n, self.n_features))

    def test_fit(self, params='stmwc', n_iter=5, verbose=False, **kwargs):
        h = hmm.GMMHMM(self.n_components)
        h.startprob = self.startprob
        h.transmat = hmm.normalize(
            self.transmat + np.diag(self.prng.rand(self.n_components)), 1)
        h.gmms = self.gmms

        # Create training data by sampling from the HMM.
        train_obs = [h.rvs(n=10, random_state=self.prng) for x in xrange(10)]

        # Mess up the parameters and see if we can re-learn them.
        h.fit(train_obs, n_iter=0)
        h.transmat = hmm.normalize(self.prng.rand(self.n_components,
                                                  self.n_components), axis=1)
        h.startprob = hmm.normalize(self.prng.rand(self.n_components))

        trainll = train_hmm_and_keep_track_of_log_likelihood(
            h, train_obs, n_iter=n_iter, params=params,
            covars_prior=1.0, **kwargs)[1:]
        if not np.all(np.diff(trainll) > 0) and verbose:
            print
            print 'Test train: (%s)\n  %s\n  %s' % (params, trainll,
                                                    np.diff(trainll))
        self.assertTrue(np.all(np.diff(trainll) > -0.5))

    def test_fit_works_on_sequences_of_different_length(self):
        obs = [self.prng.rand(3, self.n_features),
               self.prng.rand(4, self.n_features),
               self.prng.rand(5, self.n_features)]

        h = hmm.GMMHMM(self.n_components, cvtype=self.cvtype)
        # This shouldn't raise
        # ValueError: setting an array element with a sequence.
        h.fit(obs)


class TestGMMHMMWithSphericalCovars(TestGMMHMM):
    cvtype = 'spherical'

    def test_fit_startprob_and_transmat(self):
        self.test_fit('st')

    def test_fit_means(self):
        self.test_fit('m')


class TestGMMHMMWithTiedCovars(TestGMMHMM):
    cvtype = 'tied'


class TestGMMHMMWithFullCovars(TestGMMHMM):
    cvtype = 'full'
