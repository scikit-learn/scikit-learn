"""
The :mod:`sklearn.preprocessing` module includes scaling, centering,
normalization, binarization methods.
"""

from ._function_transformer import FunctionTransformer

from .data import Binarizer
from .data import KernelCenterer
from .data import MinMaxScaler
from .data import MaxAbsScaler
from .data import Normalizer
from .data import RobustScaler
from .data import StandardScaler
from .data import QuantileTransformer
from .data import add_dummy_feature
from .data import binarize
from .data import normalize
from .data import scale
from .data import robust_scale
from .data import maxabs_scale
from .data import minmax_scale
from .data import quantile_transform
from .data import power_transform
from .data import PowerTransformer
from .data import PolynomialFeatures

from ._encoders import OneHotEncoder
from ._encoders import OrdinalEncoder

from .label import label_binarize
from .label import LabelBinarizer
from .label import LabelEncoder
from .label import MultiLabelBinarizer

from ._discretization import KBinsDiscretizer


__all__ = [
    'Binarizer',
    'FunctionTransformer',
    'KBinsDiscretizer',
    'KernelCenterer',
    'LabelBinarizer',
    'LabelEncoder',
    'MultiLabelBinarizer',
    'MinMaxScaler',
    'MaxAbsScaler',
    'QuantileTransformer',
    'Normalizer',
    'OneHotEncoder',
    'OrdinalEncoder',
    'PowerTransformer',
    'RobustScaler',
    'StandardScaler',
    'add_dummy_feature',
    'PolynomialFeatures',
    'binarize',
    'normalize',
    'scale',
    'robust_scale',
    'maxabs_scale',
    'minmax_scale',
    'label_binarize',
    'quantile_transform',
    'power_transform',
]
