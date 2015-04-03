import warnings

warnings.warn("lda.LDA has been moved to "
              "discriminant_analysis.LinearDiscriminantAnalysis "
              "in 0.17 and will be removed in 0.19", DeprecationWarning)

from .discriminant_analysis import LinearDiscriminantAnalysis


class LDA(LinearDiscriminantAnalysis):
    def __init__(self, solver='svd', shrinkage=None, priors=None,
                 n_components=None, store_covariance=False, tol=1e-4):
        super(LDA, self).__init__(solver=solver, shrinkage=shrinkage, priors=priors,
                                  n_components=n_components, store_covariance=store_covariance,
                                  tol=tol)
