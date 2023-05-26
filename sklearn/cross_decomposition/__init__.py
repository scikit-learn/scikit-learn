from ..externals import _lazy_loader

__getattr__, __dir__, __all__ = _lazy_loader.attach(
    __name__,
    submod_attrs={"_pls": ["PLSSVD", "CCA", "PLSCanonical", "PLSRegression"]},
)
