from ..base import OutlierResamplerMixin
from ..covariance import EllipticEnvelope
from ..svm import OneClassSVM
from ..ensemble import IsolationForest
from ..neighbors import LocalOutlierFactor


class EllipticEnvelopeResampler(EllipticEnvelope, OutlierResamplerMixin):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class OneClassSVMResampler(OneClassSVM, OutlierResamplerMixin):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class IsolationForestResampler(IsolationForest, OutlierResamplerMixin):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class LocalOutlierFactorResampler(LocalOutlierFactor, OutlierResamplerMixin):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
