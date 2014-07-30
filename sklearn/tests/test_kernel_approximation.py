import timeit
import numpy as np
from scipy.sparse import csr_matrix
import datetime
import fht as fht
from scipy.fftpack import dct


import fht
from scipy.fftpack import dct
from scipy.linalg import hadamard

from sklearn.utils.testing import assert_array_equal, assert_equal, assert_greater
from sklearn.utils.testing import assert_not_equal, assert_almost_equal
from sklearn.utils.testing import assert_array_almost_equal, assert_raises


from sklearn.utils import check_random_state

from sklearn.metrics.pairwise import kernel_metrics
from sklearn.kernel_approximation import RBFSampler
from sklearn.kernel_approximation import AdditiveChi2Sampler
from sklearn.kernel_approximation import SkewedChi2Sampler
from sklearn.kernel_approximation import Nystroem
from sklearn.kernel_approximation import Fastfood
from sklearn.metrics.pairwise import polynomial_kernel, rbf_kernel


# generate data
rng = np.random.RandomState(0)
X = rng.random_sample(size=(300, 50))
Y = rng.random_sample(size=(300, 50))
X /= X.sum(axis=1)[:, np.newaxis]
Y /= Y.sum(axis=1)[:, np.newaxis]


def test_additive_chi2_sampler():
    """test that AdditiveChi2Sampler approximates kernel on random data"""

    # compute exact kernel
    # appreviations for easier formular
    X_ = X[:, np.newaxis, :]
    Y_ = Y[np.newaxis, :, :]

    large_kernel = 2 * X_ * Y_ / (X_ + Y_)

    # reduce to n_samples_x x n_samples_y by summing over features
    kernel = (large_kernel.sum(axis=2))

    # approximate kernel mapping
    transform = AdditiveChi2Sampler(sample_steps=3)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y)

    kernel_approx = np.dot(X_trans, Y_trans.T)

    assert_array_almost_equal(kernel, kernel_approx, 1)

    X_sp_trans = transform.fit_transform(csr_matrix(X))
    Y_sp_trans = transform.transform(csr_matrix(Y))

    assert_array_equal(X_trans, X_sp_trans.A)
    assert_array_equal(Y_trans, Y_sp_trans.A)

    # test error is raised on negative input
    Y_neg = Y.copy()
    Y_neg[0, 0] = -1
    assert_raises(ValueError, transform.transform, Y_neg)

    # test error on invalid sample_steps
    transform = AdditiveChi2Sampler(sample_steps=4)
    assert_raises(ValueError, transform.fit, X)

    # test that the sample interval is set correctly
    sample_steps_available = [1, 2, 3]
    for sample_steps in sample_steps_available:

        # test that the sample_interval is initialized correctly
        transform = AdditiveChi2Sampler(sample_steps=sample_steps)
        assert_equal(transform.sample_interval, None)

        # test that the sample_interval is changed in the fit method
        transform.fit(X)
        assert_not_equal(transform.sample_interval_, None)

    # test that the sample_interval is set correctly
    sample_interval = 0.3
    transform = AdditiveChi2Sampler(sample_steps=4,
                                    sample_interval=sample_interval)
    assert_equal(transform.sample_interval, sample_interval)
    transform.fit(X)
    assert_equal(transform.sample_interval_, sample_interval)


def test_skewed_chi2_sampler():
    """test that RBFSampler approximates kernel on random data"""

    # compute exact kernel
    c = 0.03
    # appreviations for easier formular
    X_c = (X + c)[:, np.newaxis, :]
    Y_c = (Y + c)[np.newaxis, :, :]

    # we do it in log-space in the hope that it's more stable
    # this array is n_samples_x x n_samples_y big x n_features
    log_kernel = ((np.log(X_c) / 2.) + (np.log(Y_c) / 2.) + np.log(2.) -
                  np.log(X_c + Y_c))
    # reduce to n_samples_x x n_samples_y by summing over features in log-space
    kernel = np.exp(log_kernel.sum(axis=2))

    # approximate kernel mapping
    transform = SkewedChi2Sampler(skewedness=c, n_components=1000,
                                  random_state=42)
    X_trans = transform.fit_transform(X)
    Y_trans = transform.transform(Y)

    kernel_approx = np.dot(X_trans, Y_trans.T)
    assert_array_almost_equal(kernel, kernel_approx, 1)

    # test error is raised on negative input
    Y_neg = Y.copy()
    Y_neg[0, 0] = -1
    assert_raises(ValueError, transform.transform, Y_neg)


def test_rbf_sampler():
    """test that RBFSampler approximates kernel on random data"""
    # compute exact kernel
    gamma = 10.
    kernel = rbf_kernel(X, Y, gamma=gamma)

    # approximate kernel mapping
    rbf_transform = RBFSampler(gamma=gamma, n_components=1000, random_state=42)
    X_trans = rbf_transform.fit_transform(X)
    Y_trans = rbf_transform.transform(Y)
    kernel_approx = np.dot(X_trans, Y_trans.T)


    assert_array_almost_equal(kernel, kernel_approx, 1)


def test_input_validation():
    """Regression test: kernel approx. transformers should work on lists

    No assertions; the old versions would simply crash
    """
    X = [[1, 2], [3, 4], [5, 6]]
    AdditiveChi2Sampler().fit(X).transform(X)
    SkewedChi2Sampler().fit(X).transform(X)
    RBFSampler().fit(X).transform(X)

    X = csr_matrix(X)
    RBFSampler().fit(X).transform(X)


def test_nystroem_approximation():
    # some basic tests
    rnd = np.random.RandomState(0)
    X = rnd.uniform(size=(10, 4))

    # With n_components = n_samples this is exact
    X_transformed = Nystroem(n_components=X.shape[0]).fit_transform(X)
    K = rbf_kernel(X)
    assert_array_almost_equal(np.dot(X_transformed, X_transformed.T), K)

    trans = Nystroem(n_components=2, random_state=rnd)
    X_transformed = trans.fit(X).transform(X)
    assert_equal(X_transformed.shape, (X.shape[0], 2))

    # test callable kernel
    linear_kernel = lambda X, Y: np.dot(X, Y.T)
    trans = Nystroem(n_components=2, kernel=linear_kernel, random_state=rnd)
    X_transformed = trans.fit(X).transform(X)
    assert_equal(X_transformed.shape, (X.shape[0], 2))

    # test that available kernels fit and transform
    kernels_available = kernel_metrics()
    for kern in kernels_available:
        trans = Nystroem(n_components=2, kernel=kern, random_state=rnd)
        X_transformed = trans.fit(X).transform(X)
        assert_equal(X_transformed.shape, (X.shape[0], 2))


def test_nystroem_poly_kernel_params():
    """Non-regression: Nystroem should pass other parameters beside gamma."""
    rnd = np.random.RandomState(37)
    X = rnd.uniform(size=(10, 4))

    K = polynomial_kernel(X, degree=3.1, coef0=.1)
    nystroem = Nystroem(kernel="polynomial", n_components=X.shape[0],
                        degree=3.1, coef0=.1)
    X_transformed = nystroem.fit_transform(X)
    assert_array_almost_equal(np.dot(X_transformed, X_transformed.T), K)


def test_nystroem_callable():
    """Test Nystroem on a callable."""
    rnd = np.random.RandomState(42)
    n_samples = 10
    X = rnd.uniform(size=(n_samples, 4))

    def logging_histogram_kernel(x, y, log):
        """Histogram kernel that writes to a log."""
        log.append(1)
        return np.minimum(x, y).sum()

    kernel_log = []
    X = list(X)     # test input validation
    Nystroem(kernel=logging_histogram_kernel,
             n_components=(n_samples - 1),
             kernel_params={'log': kernel_log}).fit(X)
    assert_equal(len(kernel_log), n_samples * (n_samples - 1) / 2)


def test_enforce_dimensionality_constraint():

    for message, input, expected in [
        ('test n is scaled to be a multiple of d', (16, 20), (16, 32, 2)),
        ('test n equals d', (16, 16), (16, 16, 1)),
        ('test n becomes power of two', (3, 16), (4, 16, 4)),
        ('test all', (7, 12), (8, 16, 2)),
            ]:
        d, n = input
        output = Fastfood.enforce_dimensionality_constraints(d, n)
        yield assert_equal, expected, output, message


def assert_all_equal(array_, value):
    for i in array_:
        assert_almost_equal(i, value)


def test_rows_of_gaussian_iid_have_same_length():
    fastfood = Fastfood(1, 1, 42)
    for d in (2 ** x for x in xrange(1, 8)):
        fastfood.d = d

        B, G, P, S = fastfood.create_vectors()
        HGPHB = fastfood.create_gaussian_iid_matrix(B, G, P)

        norm_of_all_rows = np.linalg.norm(HGPHB, axis=1)
        n = norm_of_all_rows[0]
        yield assert_all_equal, norm_of_all_rows, n


def test_hadamard_equivalence():
    # Ensure that the fht along axes 0 is equivalent an explicit hadamard
    # transform for random diagonal gaussian matrices
    for i, d in ((x, 2 ** x) for x in xrange(1, 8)):
        vector = np.random.normal(size=d)

        g = np.diag(np.random.normal(size=d))
        h = hadamard(d)
        normalization = (1 / np.power(2, i / 2.0))
        # explicit
        one = normalization * np.dot(h, g)
        # fht
        two = fht.fht(g, axes=0)

        # fht with vector without matrix
        a = np.dot(g,vector)
        print "a",a,g, a.shape,vector.shape
        three = fht.fht(a, axes=0)

        # fht with vector with matrix
        four = np.dot(fht.fht(g, axes=0),vector)
        
        # fht without diagonal matrix
        a = np.diag(g)*vector
        print "a",a,g, a.shape,vector.shape
        five = fht.fht(a, axes=0)
        
        yield assert_array_almost_equal, one, two
        yield assert_array_almost_equal, three, four
        yield assert_array_almost_equal, four,five


def test_scalability_to_one():
    """ Test equation 11. """
    fastfood = Fastfood(1, 1)
    d = 128
    fastfood.d = d

    B, G, P, S = fastfood.create_vectors()
    HGPHB = fastfood.create_gaussian_iid_matrix(B, G, P)

    # G is a vector need to make a diagonal matrix for the test
    G = np.diag(G)

    # The norm checks only work with an UNNORMALIZED hadamard transform
    norm_by_transpose = np.sqrt(np.diagonal(np.dot(HGPHB, HGPHB.T)))[0]
    norm_by_ax0 = np.linalg.norm(HGPHB, axis=0)[0]
    norm_by_ax1 = np.linalg.norm(HGPHB, axis=1)[0]
    norm_by_eq11 = np.sqrt(((np.linalg.norm(G) ** 2)) * d)

    assert_almost_equal(norm_by_transpose, norm_by_ax0)
    assert_almost_equal(norm_by_transpose, norm_by_ax1)
    assert_almost_equal(norm_by_transpose, norm_by_eq11)

    # This is the test for the rescale which doesn't work yet
    scale = np.power(np.linalg.norm(np.diag(G)), -1) * np.power(d, -0.5)
    print 'norm: ', norm_by_transpose
    print 'scaled: ', (norm_by_transpose) * scale
    print 'scaled matrix:', np.linalg.norm(HGPHB * scale, axis=1)


def test_V_is_gaussian():
    """test that V is a gaussian N(0,sigma^-2I_d)"""
    fastfood = Fastfood(1, 500)
    d = 512
    fastfood.d = d

    B, G, P, S = fastfood.create_vectors()
    HGPHB = fastfood.create_gaussian_iid_matrix(B, G, P)

    V = fastfood.create_approximation_matrix(S, HGPHB)

    means = np.mean(V, axis=1)
    deviations = np.std(V, axis=1)
    print 'means of rows in V', means[0:9]
    print 'deviations of rows in V', deviations[0:9]

    assert_array_almost_equal(means, np.zeros(means.shape[0]), decimal=1)
    assert_array_almost_equal(deviations, np.power(fastfood.sigma, -2)*np.ones(means.shape[0]), decimal=1)

########     Performance Analysis    #################

def test_compare_performance_of_scaling_matrix_generation():
    """comparison between scaling generation methods"""
    fastfood = Fastfood(1, 500)
    d = 2**12
    fastfood.d = d

    B, G, P, S = fastfood.create_vectors()

    start = datetime.datetime.utcnow()
    fastfood.scaling_vector_vectorized(d, G)
    end = datetime.datetime.utcnow()
    spent_time_vectorized = end - start

    start = datetime.datetime.utcnow()
    fastfood.scaling_vector_loop(d, G)
    end = datetime.datetime.utcnow()
    spent_time_loop = end - start

    start = datetime.datetime.utcnow()
    fastfood.scaling_vector_chi(d, G)
    end = datetime.datetime.utcnow()
    spent_time_chi = end - start

    print "Timimg vectorized: ", spent_time_vectorized, "Timing loop: ", spent_time_loop,"Timing chi: ", spent_time_chi

    assert_greater(spent_time_loop, spent_time_vectorized)
    assert_greater(spent_time_vectorized, spent_time_chi)


def test_fastfood():
    """test that Fastfood approximates kernel on random data"""
    # compute exact kernel
    gamma = 10.
    kernel = rbf_kernel(X, Y, gamma=gamma)

    sigma = np.sqrt(1 / (2 * gamma))

    # approximate kernel mapping
    ff_transform = Fastfood(sigma, n_components=1000, random_state=42)

    X_trans = ff_transform.fit_transform(X)
    Y_trans = ff_transform.transform(Y)
    #print X_trans, Y_trans
    kernel_approx = np.dot(X_trans, Y_trans.T)

    print 'approximation:', kernel_approx[:5, :5]
    print 'true kernel:', kernel[:5, :5]
    assert False
    assert_array_almost_equal(kernel, kernel_approx, decimal=1)

def test_fastfood_fast():
    """test that Fastfood fast approximates kernel on random data"""
    # compute exact kernel
    gamma = 10.
    kernel = rbf_kernel(X, Y, gamma=gamma)
    
    sigma = np.sqrt(1 / (2 * gamma))
    
    # approximate kernel mapping
    ff_transform = Fastfood(sigma, n_components=1000, random_state=42)
    
    pars = ff_transform.fit(X)
    X_trans = pars.transform_fast(X)
    print X_trans.shape
    Y_trans = ff_transform.transform_fast(Y)
    print Y_trans.shape


#print X_trans, Y_trans
    kernel_approx = np.dot(X_trans, Y_trans.T)
    
    print 'approximation:', kernel_approx[:5, :5]
    print 'true kernel:', kernel[:5, :5]
    assert False
    assert_array_almost_equal(kernel, kernel_approx, decimal=1)



def test_fastfood_fast_vectorized():
    """test that Fastfood fast approximates kernel on random data"""
    # compute exact kernel
    gamma = 10.
    kernel = rbf_kernel(X, Y, gamma=gamma)

    sigma = np.sqrt(1 / (2 * gamma))

    # approximate kernel mapping
    ff_transform = Fastfood(sigma, n_components=1000, random_state=42)

    pars = ff_transform.fit_vectorized(X)
    X_trans = pars.transform_fast_vectorized(X)
    print X_trans.shape
    Y_trans = ff_transform.transform_fast_vectorized(Y)
    print Y_trans.shape


    #print X_trans, Y_trans
    kernel_approx = np.dot(X_trans, Y_trans.T)

    print 'approximation:', kernel_approx[:5, :5]
    print 'true kernel:', kernel[:5, :5]
    assert False
    assert_array_almost_equal(kernel, kernel_approx, decimal=1)

def test_fastfood_mem_or_accuracy():
    """compares the performance of Fastfood and RKS"""
    #generate data
    X = rng.random_sample(size=(10000, 4000))
    X /= X.sum(axis=1)[:, np.newaxis]

    # calculate feature maps
    gamma = 10.
    sigma = np.sqrt(1 / (2 * gamma))
    number_of_features_to_generate = 1000



    fastfood_start = datetime.datetime.utcnow()
    # Fastfood: approximate kernel mapping
    rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, tradeoff_less_mem_or_higher_accuracy='accuracy', random_state=42)
    _ = rbf_transform.fit_transform(X)
    fastfood_end = datetime.datetime.utcnow()
    fastfood_spent_time =fastfood_end- fastfood_start
    print "Timimg fastfood accuracy: \t\t", fastfood_spent_time


    fastfood_mem_start = datetime.datetime.utcnow()
    # Fastfood: approximate kernel mapping
    rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, tradeoff_less_mem_or_higher_accuracy='mem', random_state=42)
    _ = rbf_transform.fit_transform(X)
    fastfood_mem_end = datetime.datetime.utcnow()
    fastfood_mem_spent_time = fastfood_mem_end- fastfood_mem_start
    print "Timimg fastfood memory: \t\t", fastfood_mem_spent_time

    assert_greater(fastfood_spent_time, fastfood_mem_spent_time)

def test_fastfood_performance_comparison_between_methods():
    """compares the performance of Fastfood and RKS"""
    #generate data
    X = rng.random_sample(size=(5000, 2000))
    Y = rng.random_sample(size=(5000, 2000))
    X /= X.sum(axis=1)[:, np.newaxis]
    Y /= Y.sum(axis=1)[:, np.newaxis]

    # calculate feature maps
    gamma = 10.
    sigma = np.sqrt(1 / (2 * gamma))
    number_of_features_to_generate = 1024


    exact_start = datetime.datetime.utcnow()
    # original rbf kernel method: 
    rbf_kernel(X, X, gamma=gamma)
    rbf_kernel(X, Y, gamma=gamma)
    exact_end = datetime.datetime.utcnow()
    exact_spent_time = exact_end- exact_start
    print "Timimg exact rbf: \t\t", exact_spent_time


    # fastfood_start = datetime.datetime.utcnow()
    # # Fastfood: approximate kernel mapping
    # rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, random_state=42)
    # _ = rbf_transform.fit_transform(X)
    # _ = rbf_transform.transform(Y)
    # fastfood_end = datetime.datetime.utcnow()
    # fastfood_spent_time =fastfood_end- fastfood_start
    # #print X_trans, Y_trans
    # #kernel_approx = np.dot(X_trans_fastfood, Y_trans_fastfood.T)
    # print "Timimg fastfood: \t\t", fastfood_spent_time
    # 
    # 
    # fastfood_fast_one_step_start = datetime.datetime.utcnow()
    # # Fastfood: approximate kernel mapping
    # rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, random_state=42)
    # #X_trans_fastfood_fast_one_step = rbf_transform.fit(X).transform_fast_one_step(X)
    # #Y_trans_fastfood_fast_one_step = rbf_transform.transform_fast_one_step(Y)
    # fastfood_fast_one_step_end = datetime.datetime.utcnow()
    # fastfood_fast_one_step_spent_time =fastfood_fast_one_step_end- fastfood_fast_one_step_start
    # #print X_trans, Y_trans
    # #kernel_approx = np.dot(X_trans_fastfood, Y_trans_fastfood.T)
    # print "Timimg fastfood fast one step: \t" ,fastfood_fast_one_step_spent_time
    # 
    # 
    # fastfood_fast_start = datetime.datetime.utcnow()
    # # Fastfood: approximate kernel mapping
    # rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, random_state=42)
    # _ = rbf_transform.fit(X).transform_fast(X)
    # _ = rbf_transform.transform_fast(Y)
    # fastfood_fast_end = datetime.datetime.utcnow()
    # fastfood_fast_spent_time =fastfood_fast_end- fastfood_fast_start
    # #print X_trans, Y_trans
    # #kernel_approx = np.dot(X_trans_fastfood, Y_trans_fastfood.T)
    # print "Timimg fastfood fast: \t\t", fastfood_fast_spent_time

    rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, 
                             tradeoff_less_mem_or_higher_accuracy='mem', random_state=42)
    _ = rbf_transform.fit_vectorized(X)
    fastfood_fast_vec_start = datetime.datetime.utcnow()
    # Fastfood: approximate kernel mapping
    rbf_transform.transform_fast_vectorized(X)
    _ = rbf_transform.transform_fast_vectorized(Y)
    fastfood_fast_vec_end = datetime.datetime.utcnow()
    fastfood_fast_vec_spent_time =fastfood_fast_vec_end- fastfood_fast_vec_start
    #print X_trans, Y_trans
    #kernel_approx = np.dot(X_trans_fastfood, Y_trans_fastfood.T)
    print "Timimg fastfood fast vectorized: \t\t", fastfood_fast_vec_spent_time


    rks_rbf_transform = RBFSampler(gamma=gamma, n_components=number_of_features_to_generate, random_state=42)
    _ = rks_rbf_transform.fit(X)
    rks_start = datetime.datetime.utcnow()
    # Random Kitchens Sinks: approximate kernel mapping
    _ = rks_rbf_transform.transform(X)
    _ = rks_rbf_transform.transform(Y)
    rks_end = datetime.datetime.utcnow()
    rks_spent_time =rks_end- rks_start
    print "Timimg rks: \t\t\t", rks_spent_time

    assert_greater(rks_spent_time, fastfood_spent_time)
    #kernel_approx = np.dot(X_trans_rks, Y_trans_rks.T)

    assert_greater(rks_spent_time, fastfood_spent_time)
    #kernel_approx = np.dot(X_trans_rks, Y_trans_rks.T)

def test_fht_dct_performance():
    """test FHT and DCT"""

    # generate data
    rng = np.random.RandomState(0)
    X = rng.random_sample(size=(4096, 4096))
    X /= X.sum(axis=1)[:, np.newaxis]

    fht_start = datetime.datetime.utcnow()
    _ = fht.fht(X)
    fht_end = datetime.datetime.utcnow()
    fht_spent_time = fht_end - fht_start

    dct_start = datetime.datetime.utcnow()
    _ = dct(X)
    dct_end = datetime.datetime.utcnow()
    dct_spent_time = dct_end- dct_start

    print "Timing fht: ", fht_spent_time, dct_spent_time

    assert_greater(fht_spent_time, dct_spent_time)

# nosetests D:\playground\scikit-learn\sklearn\tests\test_kernel_approximation.py:test_digit_recognition
def test_digit_recognition():
    print __doc__

    # Author: Gael Varoquaux <gael dot varoquaux at normalesup dot org>
    # License: Simplified BSD
    
    # Standard scientific Python imports
    import pylab as pl
    
    # Import datasets, classifiers and performance metrics
    from sklearn import datasets, svm, metrics
    from sklearn.linear_model.stochastic_gradient import SGDClassifier
    
    # The digits dataset
    digits = datasets.load_digits()
    
    # The data that we are interested in is made of 8x8 images of digits,
    # let's have a look at the first 3 images, stored in the `images`
    # attribute of the dataset. If we were working from image files, we
    # could load them using pylab.imread. For these images know which
    # digit they represent: it is given in the 'target' of the dataset.
    for index, (image, label) in enumerate(zip(digits.images, digits.target)[:4]):
        pl.subplot(2, 4, index + 1)
        pl.axis('off')
        pl.imshow(image, cmap=pl.cm.gray_r, interpolation='nearest')
        pl.title('Training: %i' % label)
    
    # To apply an classifier on this data, we need to flatten the image, to
    # turn the data in a (samples, feature) matrix:
    n_samples = len(digits.images)
    data = digits.images.reshape((n_samples, -1))
    gamma = .001
    sigma = np.sqrt(1 / (2 * gamma))
    number_of_features_to_generate = 1000
    train__idx = range(n_samples / 2)
    test__idx = range(n_samples / 2,n_samples)
    
    # Create a classifier: a support vector classifier
    classifier = svm.SVC(gamma=gamma)
    sgd_classifier = SGDClassifier()
    
    # map data into featurespace
    rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, random_state=42)
    data_transformed_train = rbf_transform.fit(data[train__idx]).transform_fast(data[train__idx])
    data_transformed_test = rbf_transform.transform_fast(data[test__idx])
    
    # We learn the digits on the first half of the digits
    classifier.fit(data[train__idx], digits.target[train__idx])
    sgd_classifier.fit(data_transformed_train, digits.target[train__idx])
    
    # Now predict the value of the digit on the second half:
    expected = digits.target[test__idx]
    predicted = classifier.predict(data[test__idx])
    predicted_sgd = sgd_classifier.predict(data_transformed_test)
    
    print "Classification report for classifier %s:\n%s\n" % (
        classifier, metrics.classification_report(expected, predicted))
    print "Classification report for classifier %s:\n%s\n" % (
        sgd_classifier, metrics.classification_report(expected, predicted_sgd))
    print "Confusion matrix:\n%s" % metrics.confusion_matrix(expected, predicted)
    
    for index, (image, prediction) in enumerate(zip(digits.images[test__idx], predicted)[:4]):
        pl.subplot(2, 4, index + 5)
        pl.axis('off')
        pl.imshow(image, cmap=pl.cm.gray_r, interpolation='nearest')
        pl.title('Prediction: %i' % prediction)
    
    pl.show()

def test_compare_fast_slow_for_correctness():

    #generate data
    X = rng.random_sample(size=(2, 2))
    #Y = rng.random_sample(size=(1, 4))
    X /= X.sum(axis=1)[:, np.newaxis]
    #Y /= Y.sum(axis=1)[:, np.newaxis]

    # calculate feature maps
    gamma = 10.
    sigma = np.sqrt(1 / (2 * gamma))
    number_of_features_to_generate = 4


    exact_start = datetime.datetime.utcnow()
    # original rbf kernel method: 
    rbf_kernel(X, X, gamma=gamma)
    #rbf_kernel(X, Y, gamma=gamma)
    exact_end = datetime.datetime.utcnow()
    exact_spent_time = exact_end- exact_start
    print "Timimg exact rbf: \t\t", exact_spent_time

    # Fastfood: approximate kernel mapping
    rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, random_state=42)
    foo = rbf_transform.fit_vectorized(X)
    foo.transform_fast(X)
    #_ = rbf_transform.transform_fast(Y)
    #print X_trans, Y_trans
    #kernel_approx = np.dot(X_trans_fastfood, Y_trans_fastfood.T)

    rbf_transform = Fastfood(sigma=sigma, n_components=number_of_features_to_generate, 
                             tradeoff_less_mem_or_higher_accuracy='mem', random_state=42)
    # Fastfood: approximate kernel mapping
    foo.transform_fast_vectorized(X)
    #_ = rbf_transform.transform_fast_vectorized(Y)
    #print X_trans, Y_trans
    #kernel_approx = np.dot(X_trans_fastfood, Y_trans_fastfood.T)
