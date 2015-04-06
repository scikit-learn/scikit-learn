import warnings
warnings.warn("qda.QDA has been moved to "
              "discriminant_analysis.QuadraticDiscriminantAnalysis "
              "in 0.17 and will be removed in 0.19.", DeprecationWarning)

from .discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
