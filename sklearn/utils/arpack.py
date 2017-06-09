# Remove this module in version 0.21

from scipy.sparse.linalg import eigs as _eigs, eigsh as _eigsh, svds as _svds

from .deprecation import deprecated


@deprecated("sklearn.utils.arpack.eigs was deprecated in version 0.19 and "
            "will be removed in 0.21. Use scipy.sparse.linalg.eigs instead.")
def eigs(A, *args, **kwargs):
    return _eigs(A, *args, **kwargs)


@deprecated("sklearn.utils.arpack.eigsh was deprecated in version 0.19 and "
            "will be removed in 0.21. Use scipy.sparse.linalg.eigsh instead.")
def eigsh(A, *args, **kwargs):
    return _eigsh(A, *args, **kwargs)


@deprecated("sklearn.utils.arpack.svds was deprecated in version 0.19 and "
            "will be removed in 0.21. Use scipy.sparse.linalg.svds instead.")
def svds(A, *args, **kwargs):
    return _svds(A, *args, **kwargs)
