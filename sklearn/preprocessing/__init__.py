"""
The :mod:`sklearn.preprocessing` module includes scaling, centering,
normalization, binarization and imputation methods.
"""

from .data import Binarizer
from .data import KernelCenterer
from .data import MinMaxScaler
from .data import Normalizer
from .data import StandardScaler
from .data import add_dummy_feature
from .data import binarize
from .data import normalize
from .data import scale

from .label import LabelBinarizer
from .label import LabelEncoder
from .label import OneHotEncoder

from .imputation import Imputer

__all__ = [
    'Binarizer',
    'Imputer',
    'KernelCenterer',
    'LabelBinarizer',
    'LabelEncoder',
    'MinMaxScaler',
    'Normalizer',
    'OneHotEncoder',
    'StandardScaler',
    'add_dummy_feature',
    'binarize',
    'normalize',
    'scale',
]