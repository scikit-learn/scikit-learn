def test():

    import numpy as np

    from sklearn.base import clone
    from sklearn.utils.testing import set_random_state

    def check_transformer_dtypes_casting(transformer, X, y):
        """Check that a trasformer preserves 64bit / 32bit
        dtypes

        Parameters
        ----------
        transformer : estimator
          a transformer instance

        X : array, shape (n_samples, n_features)
          Training vector, where n_samples is the number of samples
          and n_features is the number of features.
        y : array, shape (n_samples)
          Target vector relative to X.
        """
        for dtype_in, dtype_out in [(np.float32, np.float32),
                                    (np.float64, np.float64),
                                    (np.int, np.float64)]:
            X_cast = X.copy().astype(dtype_in)

            transformer = clone(transformer)
            set_random_state(transformer)

            if hasattr(transformer, 'fit_transform'):
                X_trans = transformer.fit_transform(X_cast, y)
            elif hasattr(transformer, 'fit_transform'):
                transformer.fit(X_cast, y)
                X_trans = transformer.transform(X_cast)

            # FIXME: should we check that the dtype of some attributes are the
            # same than dtype.
            assert X_trans.dtype == dtype_out, \
                ('transform dtype: {} - original dtype: {}'
                 .format(X_trans.dtype, X_cast.dtype))

if __name__ == "__main":

    # from sklearn.decomposition import factor_analysis
    #
    # X = np.random.RandomState(0).randn(1000, 100)
    #
    # estimator = factor_analysis.FactorAnalysis()
    #
    # check_transformer_dtypes_casting(estimator, X, None)
    #
    # # from sklearn.decomposition import dict_learning
    # #
    # # X = np.random.RandomState(0).randn(1000, 100)
    # #
    # # estimator = dict_learning.DictionaryLearning()
    # #
    # # check_transformer_dtypes_casting(estimator, X, None)
    #
    # # from sklearn.neural_network import rbm
    # #
    # # X = np.random.RandomState(0).randn(1000, 100)
    # #
    # # estimator = rbm.BernoulliRBM()
    # #
    # # check_transformer_dtypes_casting(estimator, X, None)
    #
    from sklearn.cross_decomposition import pls_

    X = np.random.RandomState(0).randn(1000, 100)

    estimator = pls_.PLSCanonical()

    check_transformer_dtypes_casting(estimator, X, None)
