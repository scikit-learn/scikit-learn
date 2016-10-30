import warnings
from .discriminant_analysis import LinearDiscriminantAnalysis as _LDA

warnings.warn("lda.LDA has been moved to "
              "discriminant_analysis.LinearDiscriminantAnalysis "
              "in 0.17 and will be removed in 0.19", DeprecationWarning)


class LDA(_LDA):
    """
    Alias for
    :class:`sklearn.discriminant_analysis.LinearDiscriminantAnalysis`.

    .. deprecated:: 0.17
        This class will be removed in 0.19.
        Use
        :class:`sklearn.discriminant_analysis.LinearDiscriminantAnalysis`
        instead.
    """
    pass
