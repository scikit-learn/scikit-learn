import pytest
import numpy as np
from sklearn.resample import EllipticEnvelopeResampler, OneClassSVMResampler, LocalOutlierFactorResampler, IsolationForestResampler
from sklearn.covariance import EllipticEnvelope
from sklearn.svm import OneClassSVM
from sklearn.ensemble import IsolationForest
from sklearn.neighbors import LocalOutlierFactor
from sklearn.datasets import make_blobs

X, y = make_blobs(random_state=0)

@pytest.mark.parametrize("detector, superclass",
                         [(EllipticEnvelopeResampler(), EllipticEnvelope()),
                          (OneClassSVMResampler(), OneClassSVM()),
                          (LocalOutlierFactorResampler(), LocalOutlierFactor()),
                          (IsolationForestResampler(), IsolationForest())])
def test_basic(detector, superclass):
    outliers = superclass.fit_predict(X, y) == -1
    n_outliers = np.sum(outliers)
    assert n_outliers > 0 # we have some outliers in the dataset

    X_new, y_new, props_new = detector.fit_resample(X, y)

    assert X_new.shape[0] == X.shape[0] - n_outliers
    assert y_new.shape[0] == y.shape[0] - n_outliers
