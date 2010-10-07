"""Feature slection module for python"""

from .univariate_selection import f_classif
from .univariate_selection import f_regression
from .univariate_selection import SelectPercentile
from .univariate_selection import SelectKBest
from .univariate_selection import SelectFpr
from .univariate_selection import SelectFdr
from .univariate_selection import SelectFwe
from .univariate_selection import GenericUnivariateSelect

from .rfe import RFE
from .rfe import RFECV

