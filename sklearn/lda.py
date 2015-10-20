import warnings
warnings.warn("lda.LDA has been moved to "
              "discriminant_analysis.LinearDiscriminantAnalysis "
              "in 0.17 and will be removed in 0.19", DeprecationWarning)

from .discriminant_analysis import LinearDiscriminantAnalysis as LDA
